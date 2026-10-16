%% A class for computing (and caching) the optimal route of each agent for a given set of tasks
%
% Used by primal_recovery_MRTA. For agent i and a set S of pending tasks (rows of the coupling
% constraints (3)), cost(i, S) solves the MILP problem of the agent with the task variables theta
% fixed to S, i.e., it computes the optimal route of the agent for performing exactly the tasks S.
% If use_dp is true, the route is first computed with Held-Karp (see route_dp), and the MILP problem
% is only solved if that route cannot be shown to be optimal.
%
% The computations are grouped in rounds (see end_round), in which the agents work in parallel. Each
% agent does its own computations of a round one after the other, so the time of a round is the
% maximum over the agents of the sum of their solver times (Section VI-C).
%

classdef RouteEvaluator < handle

    properties
        mrta % Instance of MRTA (with its agents built)
        idx_loc % m_z x m matrix: local task index (idx_Evt) of pending task j for agent i (0 if it cannot do it)
        solver_args % Cell array with name/value options for MIQPholder.solve
        cache % 1 x m cell of containers.Map: key (task set) -> struct('V', cost, 'z', solution, 'path', route)
        c_a % 1 x m cell with the cost vector (8) of each agent
        num_solves = 0 % Number of MILP problems solved
        run_time_seq = 0 % Sum of the solver times
        run_time_wall = 0 % Sum of the wall times of the calls to the solver
        run_time_par = 0 % Sum over the rounds of the time of each round
        round_times = zeros(1, 0) % 1 x m: solver time of each agent in the current round
        D % Matrix of shortest-path distances between locations
        task_v % Location of each pending task
        task_rho % Cost rho_{t_j} of each pending task
        lb_cache % 1 x m cell of containers.Map: key (set of locations) -> lower bound of the travel cost
        max_lb_locs = 12 % Maximum number of task locations for the lower bounds (Held-Karp)
        solver % AgentSolver used to solve the MILP problems of the agents
        graph = [] % RouteGraph of the location graph (used by route_dp)
        use_dp = true % Try route_dp before solving the MILP problems (only if graph is set)
        max_dp_locs = 12 % Maximum number of task locations of route_dp
        dp_tol = 1e-4 % Relative tolerance of route_dp (default relative gap tolerance of the MILP solvers)
        num_dp = 0 % Number of routes computed by route_dp
    end

    methods

        function obj = RouteEvaluator(mrta, idx_A, solver_args, agent_solver)
            % RouteEvaluator(mrta, idx_A, solver_args, [agent_solver])
            %   @idx_A: structure with fields .v and .t (location and task-type of each pending task)
            %   @solver_args: cell array with name/value options for MIQPholder.solve
            %   @agent_solver: AgentSolver used to solve the MILP problems (a sequential one if not given)
            if nargin < 4 || isempty(agent_solver)
                agent_solver = AgentSolver(mrta, solver_args, false);
            end
            obj.solver = agent_solver;
            obj.mrta = mrta;
            obj.solver_args = solver_args;
            m = mrta.m;
            m_z = length(idx_A.v);
            obj.idx_loc = zeros(m_z, m);
            obj.round_times = zeros(1, m);
            obj.cache = cell(1, m);
            obj.lb_cache = cell(1, m);
            obj.c_a = cell(1, m);
            obj.D = distances(mrta.Gv);
            obj.task_v = idx_A.v(:);
            obj.task_rho = reshape(mrta.param.cost_task(idx_A.t(:)), [], 1);
            for i = 1:m
                ag = mrta.agent{i};
                if m_z > 0
                    obj.idx_loc(:, i) = ag.idx_Evt(sub2ind(size(ag.Evt), idx_A.v(:), idx_A.t(:)));
                end
                obj.cache{i} = containers.Map('KeyType', 'char', 'ValueType', 'any');
                obj.lb_cache{i} = containers.Map('KeyType', 'char', 'ValueType', 'double');
                obj.c_a{i} = ag.build_c();
            end
        end

        function flag = can_do(self, j, i)
            % Returns true if agent @i can perform pending task @j
            flag = self.idx_loc(j, i) > 0;
        end

        function [V, z] = cost(self, i, S)
            % Returns the cost @V (Inf if infeasible) and solution @z of the optimal route of agent @i
            % for performing exactly the pending tasks @S
            key = RouteEvaluator.set_key(S);
            if ~isKey(self.cache{i}, key)
                self.evaluate(i, {S});
            end
            entry = self.cache{i}(key);
            V = entry.V;
            z = entry.z;
        end

        function evaluate(self, agents, sets, ws)
            % Computes the optimal routes of the agents agents(q) for the sets of tasks sets{q} that are
            % not cached yet (in parallel if the AgentSolver is parallel). Use cost() to get the results.
            % @ws: (optional) cell array with a feasible solution ws{q} of each problem to warm-start the
            %      solver (none if empty)
            if nargin < 4 || isempty(ws); ws = cell(1, length(agents)); end
            todo = false(1, length(agents));
            keys_q = cell(1, length(agents));
            batch = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            for q = 1:length(agents)
                keys_q{q} = RouteEvaluator.set_key(sets{q});
                batch_key = sprintf('%d:%s', agents(q), keys_q{q});
                if ~isKey(self.cache{agents(q)}, keys_q{q}) && ~isKey(batch, batch_key)
                    todo(q) = true;
                    batch(batch_key) = true;
                end
            end
            if ~any(todo); return; end

            % Try first without solving the MILP problems
            if self.use_dp && ~isempty(self.graph)
                for q = find(todo)
                    i = agents(q);
                    t_start = tic;
                    [ok, z_q, V_q, path_q] = self.route_dp(i, sets{q});
                    if ok
                        self.cache{i}(keys_q{q}) = struct('V', V_q, 'z', z_q, 'path', path_q);
                        todo(q) = false;
                        self.num_dp = self.num_dp + 1;
                    end
                    t_q = toc(t_start); % Done by the agent, so it counts in the current round
                    self.run_time_seq = self.run_time_seq + t_q;
                    self.run_time_wall = self.run_time_wall + t_q;
                    self.round_times(i) = self.round_times(i) + t_q;
                end
                if ~any(todo); return; end
            end
            agents = agents(todo);
            sets = sets(todo);
            keys_q = keys_q(todo);
            ws = ws(todo);

            n = length(agents);
            c = cell(1, n);
            lbz = cell(1, n);
            ubz = cell(1, n);
            for q = 1:n
                i = agents(q);
                c{q} = self.c_a{i};
                [lbz{q}, ubz{q}] = self.mrta.agent{i}.fixed_task_bounds(self.idx_loc(sets{q}, i));
            end
            [z, ~, info, wall_time] = self.solver.solve(agents, c, lbz, ubz, ws);

            for q = 1:n
                i = agents(q);
                if info{q}.feasible
                    z_q = self.mrta.agent{i}.round_z(z{q});
                    V_q = self.c_a{i}'*z_q;
                    path_q = self.mrta.agent{i}.get_path(z_q);
                else
                    z_q = [];
                    V_q = inf;
                    path_q = [];
                end
                self.cache{i}(keys_q{q}) = struct('V', V_q, 'z', z_q, 'path', path_q);
                self.add_time(i, info{q});
            end
            self.run_time_wall = self.run_time_wall + wall_time;
        end

        function V = lower_bound(self, i, S)
            % Returns a lower bound of cost(i, S): the cost if it is cached; otherwise, the length of the
            % shortest walk from the location of the agent that visits the locations of the tasks S
            % (Held-Karp on the shortest-path distances), plus the stopping cost and the task costs.
            % Returns -Inf if the tasks are at more than max_lb_locs different locations.
            S = sort(S(:))';
            key = RouteEvaluator.set_key(S);
            if isKey(self.cache{i}, key)
                entry = self.cache{i}(key);
                V = entry.V;
                return
            end
            ag = self.mrta.agent{i};
            locs = setdiff(unique(self.task_v(S)), ag.loc);
            if length(locs) > self.max_lb_locs
                V = -inf;
                return
            end
            loc_key = ['l', sprintf('%d,', locs)];
            if isKey(self.lb_cache{i}, loc_key)
                travel = self.lb_cache{i}(loc_key);
            else
                travel = RouteEvaluator.open_walk_length(self.D, ag.loc, locs);
                self.lb_cache{i}(loc_key) = travel;
            end
            V = travel + ag.param.cost_stop + sum(self.task_rho(S));
        end

        function [ok, z, V, route] = route_dp(self, i, S)
            % Computes the route of agent @i for the pending tasks @S without solving its MILP problem
            % (Section VI-B): the shortest walk from its location that visits the locations of the tasks
            % (Held-Karp on the shortest-path distances). If the walk is a simple path and its cost is
            % within dp_tol of the lower bound
            %   LB = length of the walk + stopping cost + task costs
            %        + epsilon_tau*(sum of the shortest-path distances to the task locations),
            % the route is optimal up to the relative gap tolerance of the MILP solvers, and ok = true is
            % returned with the solution @z, its cost @V and the route. Otherwise, ok = false.
            % The lower bound assumes that the stopping cost is the same at every location.
            ok = false;
            z = [];
            V = inf;
            ag = self.mrta.agent{i};
            s = ag.loc;
            locs = reshape(setdiff(unique(self.task_v(S)), s), 1, []);
            k = length(locs);
            if k > self.max_dp_locs
                route = [];
                return
            end
            len = 0;
            route = s;
            if k > 0
                [len, order] = RouteEvaluator.open_walk_order(self.D, s, locs);
                if ~isfinite(len); return; end
                for t = 1:k
                    leg = self.graph.sp_path(route(end), locs(order(t)));
                    route = [route, leg(2:end)]; %#ok<AGROW>
                end
                if length(unique(route)) < length(route); return; end % Not a simple path
            end
            z_r = ag.route_to_z(route, self.idx_loc(S, i));
            V_r = self.c_a{i}'*z_r;
            LB = len + ag.param.cost_stop + sum(self.task_rho(S)) + ag.param.cost_times*sum(self.D(s, locs));
            if V_r - LB <= self.dp_tol*max(1, abs(V_r)) && ag.miqp.is_feasible(z_r)
                ok = true;
                z = z_r;
                V = V_r;
            end
        end

        function V = store(self, i, S, z)
            % Stores the solution @z of agent @i for the tasks @S (computed elsewhere) in the cache
            V = self.c_a{i}'*z;
            self.cache{i}(RouteEvaluator.set_key(S)) = struct('V', V, 'z', z, 'path', self.mrta.agent{i}.get_path(z));
        end

        function record_time(self, i, info)
            % Records the time of a MILP problem solved by agent @i outside of evaluate/price_agents
            self.add_time(i, info);
            self.run_time_wall = self.run_time_wall + info.wall_time;
        end

        function [V, bound, z, info] = price_agents(self, duals)
            % Pricing problems of the column generation for the dual variables @duals (one per pending
            % task): MILP problem of each agent with cost rho_{t_j} - duals(j) for each pending task j.
            % Returns the optimal values @V, lower bounds @bound of them, and the solutions @z and @info.
            % The bounds of Agent.prize_bounds are used if the agents use the valid inequalities.
            m = length(self.c_a);
            c = cell(1, m);
            lbz = cell(1, m);
            ubz = cell(1, m);
            for i = 1:m
                ag = self.mrta.agent{i};
                theta_cost = self.c_a{i}(ag.dv.n_e + ag.dv.n + (1:ag.dv.n_t));
                j = find(self.idx_loc(:, i) > 0);
                theta_cost(self.idx_loc(j, i)) = theta_cost(self.idx_loc(j, i)) - duals(j);
                c{i} = ag.task_cost_vector(theta_cost);
                if ag.tight
                    [lbz{i}, ubz{i}] = ag.prize_bounds(c{i}, ag.miqp.primal.lbz, ag.miqp.primal.ubz);
                end
            end
            [z, V, info, wall_time] = self.solver.solve(1:m, c, lbz, ubz);
            bound = zeros(1, m);
            for i = 1:m
                bound(i) = info{i}.obj_bound;
                if info{i}.feasible
                    z{i} = self.mrta.agent{i}.round_z(z{i});
                end
                self.add_time(i, info{i});
            end
            self.run_time_wall = self.run_time_wall + wall_time;
        end

        function end_round(self)
            % Ends a round of computations done in parallel by the agents
            self.run_time_par = self.run_time_par + max([self.round_times, 0]);
            self.round_times = zeros(1, length(self.c_a));
        end

        function paths = get_paths(self, i)
            % Returns a cell array with the routes of all the cached solutions of agent @i
            vals = values(self.cache{i});
            paths = cell(1, 0);
            for c = 1:length(vals)
                if ~isempty(vals{c}.path)
                    paths{end+1} = vals{c}.path; %#ok<AGROW>
                end
            end
        end

    end

    methods (Access = private)

        function add_time(self, i, info)
            % Records the solver time of a MILP problem solved by agent @i in the current round
            self.num_solves = self.num_solves + 1;
            self.run_time_seq = self.run_time_seq + info.run_time;
            self.round_times(i) = self.round_times(i) + info.run_time;
        end

    end

    methods (Static)

        function key = set_key(S)
            % Key of the set of tasks @S in the cache
            key = ['s', sprintf('%d,', sort(S(:))')];
        end

        function [len, order] = open_walk_order(D, start, locs)
            % Held-Karp: length of the shortest walk that starts at @start and visits all the locations
            % @locs (see open_walk_length), and the order in which it visits them (indices of @locs)
            k = length(locs);
            n_m = 2^k;
            Dl = D(locs, locs);
            bits = 2.^(0:k-1);
            dp = inf(n_m, k); % dp(mask+1, q): shortest walk visiting the set mask and ending at locs(q)
            par = zeros(n_m, k); % Previous location of the walk of dp(mask+1, q)
            dp(bits + 1 + (0:k-1)*n_m) = D(start, locs);
            for mask = 1:(n_m - 1)
                cur = dp(mask + 1, :);
                if all(isinf(cur)); continue; end
                [best, arg] = min(cur' + Dl, [], 1); % Shortest walk visiting mask and then each location
                out = find(bitand(mask, bits) == 0); % Locations not in mask
                idx = mask + bits(out) + 1 + (out - 1)*n_m;
                upd = best(out) < dp(idx);
                dp(idx(upd)) = best(out(upd));
                par(idx(upd)) = arg(out(upd));
            end
            [len, last] = min(dp(n_m, :));
            order = zeros(1, k);
            mask = n_m - 1;
            for t = k:-1:1
                order(t) = last;
                q = par(mask + 1, last);
                mask = mask - bits(last);
                last = q;
            end
        end

        function len = open_walk_length(D, start, locs)
            % Held-Karp: length of the shortest walk that starts at @start and visits all the locations
            % @locs, where D(u, v) is the shortest-path distance from u to v
            k = length(locs);
            if k == 0
                len = 0;
                return
            end
            Dl = D(locs, locs);
            bits = 2.^(0:k-1);
            dp = inf(2^k, k); % dp(mask+1, q): shortest walk visiting the set mask and ending at locs(q)
            dp(bits + 1 + (0:k-1)*2^k) = D(start, locs);
            for mask = 1:(2^k - 1)
                cur = dp(mask + 1, :);
                if all(isinf(cur)); continue; end
                best = min(cur' + Dl, [], 1); % Shortest walk visiting mask and then each location
                out = find(bitand(mask, bits) == 0); % Locations not in mask
                idx = mask + bits(out) + 1 + (out - 1)*2^k;
                dp(idx) = min(dp(idx), best(out));
            end
            len = min(dp(end, :));
        end

    end

end
