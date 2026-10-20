%% A class for defining the agents (robots) in the MRTA class
%
% The class contains the properties of the agent, such as the tasks it can
% perform, as well as methods that are useful for building the MRTA problem.
%
% The decision variables z_i of agent s_i (Section III of the article) are ordered as
%   z_i = [x_{i,j} (edges, n_e), x_{i,i} (stops, n), theta_{l,k} (tasks, n_t), tau_l (times, n)]
% and its MILP problem (set Z_i in (11)) is given by constraints (1), (4), (5), (6), and (7) if
% 'tight' is true.
%

classdef Agent < handle

    properties
        Gv % Location graph of the agent
        name % Name of the agent
        loc % Location of the agent
        l_with_t % Array containing the indices of locations where the agent can perform at least one task
        n_l_with_t % Number of locations where the agent can perform at least one task
        ts_at_l % Array containing the number of tasks the agent can perform at each location, |T_s(v_l)|
        t_loc % Local copy of the current pending tasks
        Eat % Array defining the tasks it can perform (row of the global Eat for this agent)
        Evt % Matrix defining the tasks it can perform in each location (takes into account if the task can be assigned to that location)
        idx_Evt % Index of each location-task pair in the part of the decision variables that encodes the tasks
        task_idx % Indices of the locations and tasks that the agent can do (the rows and cols of Evt that are 1)
        dv % Struct defining the dimensions associated with its decision variables
           % dv.idx - Integer that identifies the index in the decision variables corresponding to this agent
           % dv.n_e - Number of edges that the agent can travel through
           % dv.n - Number of locations where the agent can stop
           % dv.n_t - Number of tasks-location the agent can do
           % dv.nz - Total number of decision variables associated with the agent

        % For defining and working with its own MILP problem
        miqp % Instance of MIQPholder that defines the agent's MILP problem (13)
        param % Structure containing information like the cost of each task, etc.
        tight = false % If true, the constraints (7) are added to the MILP problem (see gen_valid_ineq_constraints)

    end

    properties (Constant)
        max_cycle_sets = 2000 % Maximum number of sets (cycles of 3 or 4 locations) used by gen_valid_ineq_constraints
    end

    methods

        % Constructor
        function obj = Agent(Gv, Evt, Eat, opt)
            arguments
                Gv digraph
                Evt (:, :) % Matrix defining the tasks that can be assigned to each location
                Eat (:, :) % Array defining the tasks it can perform (row of the global Eat for this agent)
                opt.name string = "" % Name of the agent
                opt.loc (1, 1) = 0 % Current location of the agent
                opt.t_loc = {} % Local copy of the current pending tasks
                opt.param struct = struct() % Contains parameters of the MRTA problem, such as costs
                opt.idx (1, 1) = 0 % Index where the agent starts in the global decision vector
                opt.tight (1, 1) {mustBeNumericOrLogical} = false % Add valid inequalities (see gen_valid_ineq_constraints)
            end
            obj.Gv = Gv;
            obj.tight = logical(opt.tight);
            obj.name = opt.name;
            obj.param = opt.param;
            obj.loc = opt.loc;
            obj.t_loc = opt.t_loc;
            obj.Eat = Eat;
            obj.Evt = Eat.*Evt;
            obj.dv.idx = opt.idx;
            obj.dv.n_e = numedges(obj.Gv);
            obj.dv.n = numnodes(obj.Gv);
            obj.dv.n_t = nnz(obj.Evt);
            obj.dv.nz = obj.dv.n_e + 2*obj.dv.n + obj.dv.n_t;

            % Index of each task-location pair (ordered by location, then by task)
            obj.idx_Evt = zeros(size(obj.Evt));
            [t_col, v_row] = find(obj.Evt');
            obj.idx_Evt(sub2ind(size(obj.Evt), v_row, t_col)) = 1:obj.dv.n_t;

            obj.ts_at_l = sum(obj.Evt, 2);
            obj.l_with_t = find(obj.ts_at_l > 0)';
            obj.n_l_with_t = length(obj.l_with_t);

            [obj.task_idx.v, obj.task_idx.t] = find(obj.Evt == 1);

        end

        function build_MILP(self)
            % Builds the MILP problem of the agent (without the coupling constraints (3)) in self.miqp

            % Cost function (8)
            c = self.build_c();

            % Equality constraints (1)
            [Cio, dio] = self.gen_in_out_eq_constraints();

            % Inequality constraints (4) and (5)
            [Ccd, lcd, ucd] = self.gen_can_do_ineq_constraints();
            [Cst, lst, ust] = self.gen_subtour_ineq_constraints();

            % Box constraints (6) and binary variables (9c), (9d)
            [lbz, ubz] = self.gen_box_constraints();
            [D, li, ui] = self.gen_integer_constraints();

            % Construct MILP problem
            Q = sparse(self.dv.nz, self.dv.nz); % MILP problem
            A = [Ccd; Cst];
            lb = [lcd; lst];
            ub = [ucd; ust];
            if self.tight
                [Cv, lv, uv] = self.gen_valid_ineq_constraints();
                A = [sparse(A); Cv];
                lb = [lb; lv];
                ub = [ub; uv];
            end
            idx_int = self.get_idx_integrality();
            self.miqp = MIQPholder(Q, c, A, lb, ub, Cio, dio, D, li, ui, lbz, ubz, idx_int);

        end

        function c = build_c(self)
            % Returns the vector c with the costs (8) of the decision variables of the agent
            c = zeros(self.dv.nz, 1);
            c(1:self.dv.n_e) = self.Gv.Edges.Weight; % w_{i,j}: move between locations
            c(self.dv.n_e + (1:self.dv.n)) = self.param.cost_stop; % w_{i,i}: stopping at a location
            w_tasks = (self.param.cost_task.*self.Evt)'; % rho_{t_k}: completing a task
            w_tasks = w_tasks(self.Evt' == 1); % Extract elements in the dv vector (ordered as idx_Evt)
            c(self.dv.n_e + self.dv.n + (1:self.dv.n_t)) = w_tasks;
            c(self.dv.n_e + self.dv.n + self.dv.n_t + (1:self.dv.n)) = self.param.cost_times; % epsilon_tau: time variables
        end

        function [Cio, dio] = gen_in_out_eq_constraints(self)
            % Returns matrix Cio and vector dio of the equality constraints (1):
            %   (outgoing edges + stop) - (incoming edges) = delta_l, at each location l
            Cio = zeros(self.dv.n, self.dv.nz);
            for l = 1:self.dv.n
                Cio(l, self.Gv.inedges(l)) = -1;
                Cio(l, self.Gv.outedges(l)) = 1;
                Cio(l, self.dv.n_e + l) = 1;
            end
            Cio = sparse(Cio);
            dio = self.update_dio();
        end

        function dio = update_dio(self, in_loc)
            % Updates the vector dio of gen_in_out_eq_constraints for the current agent location @in_loc.
            % If no @in_loc is provided, then the agent't self.loc variable is used.
            % If an @in_loc is provided, then the agent updates its self.loc to the new location.
            arguments
                self Agent
                in_loc (1, 1) = NaN
            end
            if isnan(in_loc)
                in_loc = self.loc;
            else
                self.loc = in_loc;
            end
            dio = zeros(self.dv.n, 1);
            dio(in_loc) = 1;
        end

        function [Ccd, lcd, ucd] = gen_can_do_ineq_constraints(self)
            % Returns matrix Ccd and vectors lcd, ucd of the inequality constraints (4),
            % divided by |T_s(v_l)|: lcd <= Ccd*z <= ucd, i.e.,
            %   (1/|T_s(v_l)|)*sum_k theta_{l,k} - sum_i x_{i,l} <= delta_l
            Ccd = zeros(self.n_l_with_t, self.dv.nz);
            k = size(self.Evt, 2);
            for v = 1:self.n_l_with_t
                l = self.l_with_t(v);
                Ccd(v, self.Gv.inedges(l)) = -1;
                for t = 1:k
                    [has_task, ~, idx_tv] = self.has_task_loc(t, l);
                    if has_task
                        Ccd(v, idx_tv) = 1/self.ts_at_l(l);
                    end
                end
            end
            Ccd = sparse(Ccd);
            [lcd, ucd] = self.update_cd_bounds();
        end

        function [lcd, ucd] = update_cd_bounds(self, in_loc)
            % Updates the lower and upper bounds lcd and ucd of gen_can_do_ineq_constraints() for the current agent location @in_loc.
            % If no @in_loc is provided, then the agent's self.loc variable is used.
            % If an @in_loc is provided, then the agent updates its self.loc to the new location.
            arguments
                self Agent
                in_loc (1, 1) = NaN
            end
            if isnan(in_loc)
                in_loc = self.loc;
            else
                self.loc = in_loc;
            end
            ucd = zeros(self.n_l_with_t, 1);
            if self.ts_at_l(in_loc) > 0
                ucd(self.l_with_t == in_loc) = 1;
            end
            lcd = -inf*ones(self.n_l_with_t, 1);
        end

        function L = big_L(self)
            % Returns L_s = sum of the edge weights w_{i,j} (used in (5) and as upper bound of tau)
            L = sum(self.Gv.Edges.Weight);
        end

        function [Cst, lst, ust] = gen_subtour_ineq_constraints(self, opt)
            % Returns matrix Cst and vectors lst, ust of the subtour elimination constraints (5):
            %   tau_j >= tau_i + w_{i,j} - L_s*(1 - x_{i,j})  <=>  L_s*x_{i,j} + tau_i - tau_j <= L_s - w_{i,j}
            % If 'with_box' is true, the box constraints 0 <= tau <= L_s are added as inequality rows.
            arguments
                self Agent
                opt.with_box (1, 1) {mustBeNumericOrLogical} = false;
            end

            L = self.big_L();
            [s_node, t_node] = findedge(self.Gv);
            e = (1:self.dv.n_e)';
            idx_t_var = self.dv.n_e + self.dv.n + self.dv.n_t;
            Cst = sparse([e; e; e], [e; idx_t_var + s_node; idx_t_var + t_node], ...
                         [L*ones(self.dv.n_e, 1); ones(self.dv.n_e, 1); -ones(self.dv.n_e, 1)], ...
                         self.dv.n_e, self.dv.nz);
            lst = -inf*ones(self.dv.n_e, 1);
            ust = L - self.Gv.Edges.Weight;

            if opt.with_box
                Cst = [Cst; sparse(1:self.dv.n, idx_t_var + (1:self.dv.n), 1, self.dv.n, self.dv.nz)];
                lst = [lst; zeros(self.dv.n, 1)];
                ust = [ust; L*ones(self.dv.n, 1)];
            end

        end

        function [lbz, ubz] = gen_box_constraints(self)
            % Returns the box constraints lbz <= z <= ubz.
            % Binary variables are in [0, 1], and the times satisfy (6), tau >= 0. The upper bound
            % tau <= L_s is not in the article, but every optimal solution satisfies it.
            lbz = zeros(self.dv.nz, 1);
            ubz = [ones(self.dv.n_e + self.dv.n + self.dv.n_t, 1); self.big_L()*ones(self.dv.n, 1)];
            if self.tight
                % The route never enters the start location, whose time is 0 (7b)
                ubz(self.Gv.inedges(self.loc)) = 0;
                ubz(self.dv.n_e + self.dv.n + self.dv.n_t + self.loc) = 0;
            end
        end

        function [Cv, lv, uv] = gen_valid_ineq_constraints(self)
            % Returns matrix Cv and vectors lv, uv of the constraints (7), lv <= Cv*z <= uv (used if
            % self.tight is true). Every route of the agent satisfies them, so they do not change the optimal
            % routes, but they make the LP relaxation much tighter. With y_v = sum_u x_{u,v} (v ~= self.loc)
            % and y_{self.loc} = 1:
            %   - y_v <= 1: a location is entered at most once (7a)
            %   - theta_{v,t} <= y_v for each task t at location v ~= self.loc (constraint (4) for each task)
            %   - x_{i,j} + x_{j,i} <= y_k for k = i, j: (7c) for the pairs of locations {i, j}
            %   - sum_{(i,j) in E(S)} x_{i,j} <= sum_{v in S, v ~= k} y_v for each k in S: (7c) for the sets S
            %     of 3 or 4 locations that form a cycle of the location graph (only k = self.loc if self.loc
            %     is in S)
            % The bounds x_{i,self.loc} = 0 and tau_{self.loc} = 0 (7b) are set in gen_box_constraints.
            % The constraints depend on self.loc, so the MILP problem must be rebuilt if the agent moves.
            n = self.dv.n;
            n_e = self.dv.n_e;
            [s_node, t_node] = findedge(self.Gv);
            s_node = s_node(:);
            t_node = t_node(:);
            loc = self.loc;
            in_e = cell(1, n); % Edges entering each location
            for v = 1:n
                in_e{v} = reshape(self.Gv.inedges(v), 1, []);
            end
            % Each row is stored as {columns, coefficients, upper bound}; -y_v is written as -1 on the edges
            % entering v (and y_loc = 1 as a constant on the right-hand side)
            rows = cell(0, 3);

            % y_v <= 1 and theta_{v,t} <= y_v
            for v = [1:(loc-1), (loc+1):n]
                if isempty(in_e{v}); continue; end
                rows(end+1, :) = {in_e{v}, ones(size(in_e{v})), 1}; %#ok<AGROW>
                for t = find(self.idx_Evt(v, :))
                    rows(end+1, :) = {[n_e + n + self.idx_Evt(v, t), in_e{v}], [1, -ones(size(in_e{v}))], 0}; %#ok<AGROW>
                end
            end

            % Pairs {i, j} connected in both directions
            pair = sparse(s_node, t_node, 1:n_e, n, n);
            for e = 1:n_e
                i = s_node(e); j = t_node(e);
                e2 = full(pair(j, i));
                if i < j && e2 > 0
                    for k = [i, j]
                        if k == loc; continue; end
                        rows(end+1, :) = {[e, e2, in_e{k}], [1, 1, -ones(size(in_e{k}))], 0}; %#ok<AGROW>
                    end
                end
            end

            % Sets S of 3 or 4 locations that form a cycle of the (undirected) location graph
            sets = Agent.short_cycles(sparse(s_node, t_node, 1, n, n), Agent.max_cycle_sets);
            for q = 1:size(sets, 1)
                S = sets(q, sets(q, :) > 0);
                in_S = false(n, 1);
                in_S(S) = true;
                arcs = find(in_S(s_node) & in_S(t_node))';
                if any(S == loc)
                    ks = loc;
                else
                    ks = S;
                end
                for k = ks
                    rest = S(S ~= k);
                    cols = [in_e{rest(rest ~= loc)}];
                    rows(end+1, :) = {[arcs, cols], [ones(size(arcs)), -ones(size(cols))], double(any(rest == loc))}; %#ok<AGROW>
                end
            end

            r = size(rows, 1);
            lv = -inf(r, 1);
            if r == 0
                Cv = sparse(0, self.dv.nz);
                uv = zeros(0, 1);
                return
            end
            len = cellfun(@numel, rows(:, 1));
            I = repelem((1:r)', len(:));
            J = [rows{:, 1}]';
            W = [rows{:, 2}]';
            Cv = sparse(I, J, W, r, self.dv.nz); % Repeated indices are added (e.g., an edge of E(S) entering v)
            uv = [rows{:, 3}]';
        end

        function [lbz, ubz] = prize_bounds(self, c, lbz, ubz)
            % Returns the bounds @lbz, @ubz with some variables fixed to 0 without changing the optimal value
            % of the MILP problem of the agent with costs @c (Section VII): the tasks with cost c >= 0, and
            % the stops at locations other than self.loc without any task of cost c < 0 (stopping there is
            % more expensive than stopping one edge earlier, since w_{i,j} > 0 and all the stopping costs
            % are equal). @lbz and @ubz are returned unchanged if some edge weight is not positive.
            if any(self.Gv.Edges.Weight <= 0); return; end
            idx_theta = self.dv.n_e + self.dv.n + (1:self.dv.n_t);
            neg = c(idx_theta) < 0;
            ubz(idx_theta(~neg)) = 0;
            lbz(idx_theta(~neg)) = 0;
            stop_ok = false(self.dv.n, 1);
            stop_ok(self.theta_loc(neg)) = true;
            stop_ok(self.loc) = true;
            ubz(self.dv.n_e + find(~stop_ok)) = 0;
            lbz(self.dv.n_e + find(~stop_ok)) = 0;
        end

        function [D, li, ui] = gen_integer_constraints(self)
            % Returns matrix D and vectors li, ui that determine the integrality constraints
            n_int = self.dv.n_e + self.dv.n + self.dv.n_t;
            D = [speye(n_int), sparse(n_int, self.dv.n)];
            li = zeros(n_int, 1);
            ui = ones(n_int, 1);
        end

        function int_idx = get_idx_integrality(self)
            % Returns the indices of the integer (binary) variables of the agent's MILP problem
            int_idx = 1:(self.dv.n_e + self.dv.n + self.dv.n_t); % The time variables tau are the last ones
        end

        function [has, idx, idx_local] = has_task_loc(self, t, v)
            % Returns true if the agent has the given task @t at location @v
            % Returns the index of the decision variables associated to the task-location pair.
            % Also returns the local index (instead of the one related to the complete dv vector).
            % Returns an empty [] if the agent does not have that task associated to that location
            if self.idx_Evt(v, t) ~= 0
                has = true;
                idx_local = self.dv.n_e + self.dv.n + self.idx_Evt(v, t);
                idx = idx_local + self.dv.idx;
            else
                has = false;
                idx_local = [];
                idx = [];
            end
        end

        function z = round_z(self, z)
            % Returns @z where the binary elements have been rounded to 0 or 1
            idx = 1:(self.dv.n_e + self.dv.n + self.dv.n_t);
            z(idx) = double(z(idx) >= 0.5);
        end

        function [path, edges] = get_path(self, z)
            % Returns the ordered list of locations @path visited by the agent (starting at self.loc)
            % and the indices @edges of the edges it travels through, according to the solution @z
            path = self.loc;
            edges = zeros(1, 0);
            cur = self.loc;
            for step = 1:self.dv.n
                out = self.Gv.outedges(cur);
                e = out(z(out) >= 0.5);
                if isempty(e); break; end
                [~, cur] = findedge(self.Gv, e(1));
                path(end+1) = cur; %#ok<AGROW>
                edges(end+1) = e(1); %#ok<AGROW>
            end
        end

        function Gam = task_times(self, z)
            % Returns the n x k matrix @Gam with the times Gamma_i(t_j) at which the agent completes each task
            % in the solution @z (Section V-B): weights w_{i,j} of the edges travelled plus the times
            % rho_{t_k} of the tasks performed so far, in the order of the path.
            % Gam(v, t) = 0 if the agent cannot do task t at location v, and Inf if it can but does not.
            Gam = self.Evt;
            Gam(Gam == 1) = inf;
            theta = z(self.dv.n_e + self.dv.n + (1:self.dv.n_t));
            [path, edges] = self.get_path(z);
            time_acc = 0;
            for p = 1:length(path)
                if p > 1
                    time_acc = time_acc + self.Gv.Edges.Weight(edges(p-1));
                end
                v = path(p);
                for t = find(self.Evt(v, :))
                    if theta(self.idx_Evt(v, t)) >= 0.5
                        time_acc = time_acc + self.param.cost_task(t);
                        Gam(v, t) = time_acc;
                    end
                end
            end
        end

        function [lbz, ubz] = fixed_task_bounds(self, tasks_on)
            % Returns the bounds of the decision variables of the agent's MILP problem with the task
            % variables theta fixed: theta = 1 for the local task indices @tasks_on (indices of idx_Evt),
            % and theta = 0 otherwise.
            % If self.tight is true (and the edge weights are positive), the stops at locations other than
            % self.loc and the locations of the tasks are also fixed to 0, since the optimal route ends at
            % the location of one of its tasks.
            lbz = self.miqp.primal.lbz;
            ubz = self.miqp.primal.ubz;
            idx = self.dv.n_e + self.dv.n + (1:self.dv.n_t);
            lbz(idx) = 0;
            ubz(idx) = 0;
            lbz(idx(tasks_on)) = 1;
            ubz(idx(tasks_on)) = 1;
            if self.tight && all(self.Gv.Edges.Weight > 0)
                stop_ok = false(self.dv.n, 1);
                stop_ok(self.theta_loc(tasks_on)) = true;
                stop_ok(self.loc) = true;
                ubz(self.dv.n_e + find(~stop_ok)) = 0;
            end
        end

        function v = theta_loc(self, q)
            % Returns the locations of the task variables theta with local indices @q (indices of idx_Evt)
            [v_list, t_list] = find(self.idx_Evt > 0);
            v_all = zeros(self.dv.n_t, 1);
            v_all(self.idx_Evt(sub2ind(size(self.idx_Evt), v_list, t_list))) = v_list;
            v = v_all(q);
        end

        function c = task_cost_vector(self, theta_cost)
            % Returns the costs (8) of the agent with the costs of the task variables theta replaced by
            % @theta_cost (vector of length dv.n_t)
            c = self.build_c();
            c(self.dv.n_e + self.dv.n + (1:self.dv.n_t)) = theta_cost;
        end

        function z = route_to_z(self, route, tasks_on)
            % Returns the solution z of the MILP problem of the agent that travels the route @route (simple
            % path starting at self.loc), stops at its last location and performs the tasks with local
            % indices @tasks_on (indices of idx_Evt). The times tau are the accumulated travel times at the
            % visited locations (0 at the other locations).
            route = reshape(route, 1, []);
            z = zeros(self.dv.nz, 1);
            if length(route) > 1
                e = findedge(self.Gv, route(1:end-1), route(2:end));
                z(e) = 1;
                z(self.dv.n_e + self.dv.n + self.dv.n_t + route(2:end)) = cumsum(self.Gv.Edges.Weight(e));
            end
            z(self.dv.n_e + route(end)) = 1;
            z(self.dv.n_e + self.dv.n + tasks_on) = 1;
        end

        function z = warm_start_z(self, z_prev, c)
            % Returns a feasible solution of the MILP problem of the agent with costs @c (and the bounds of
            % prize_bounds) built from its solution @z_prev at the previous iteration (Section VII): the
            % route of @z_prev up to its last location with a task of negative cost, performing the tasks
            % of negative cost along it. Used to warm-start the solver.
            route = self.get_path(z_prev);
            idx_theta = self.dv.n_e + self.dv.n + (1:self.dv.n_t);
            neg = c(idx_theta) < 0;
            theta_loc = self.theta_loc(1:self.dv.n_t);
            has_neg = false(self.dv.n, 1);
            has_neg(theta_loc(neg)) = true;
            last = find(has_neg(route), 1, 'last');
            if isempty(last); last = 1; end
            route = route(1:last);
            on_route = false(self.dv.n, 1);
            on_route(route) = true;
            tasks_on = find(neg(:) & on_route(theta_loc));
            z = self.route_to_z(route, tasks_on);
        end

        function text = res2tex(self, z, Gv, names)
            % Generates the text that interprets the solution @z of this agent.
            % To print human-interpretable text, the function requires the location graph @Gv
            % and the @names structure containing the names of locations and tasks
            % @z is assumed to be rounded so that binary variables are exactly 0 or 1

            text = "";
            text = strcat(text, "Agent ", self.name, "\n");
            text = strcat(text, "- Starting location: ", names.v{self.loc}, "\n");

            % Movements
            text = strcat(text, "- Movements:");
            if sum(z(1:self.dv.n_e)) == 0
                text = strcat(text, " none\n");
            else
                text = strcat(text, "\n");
                for i = 1:self.dv.n_e
                    if z(i) == 1
                        [origin, destination] = findedge(Gv, i);
                        text = strcat(text, "  > From ", names.v{origin}, " to ", names.v{destination}, "\n");
                    end
                end
            end
            text = strcat(text, "- Stops at: ");
            for i = 1:self.dv.n
                if z(self.dv.n_e + i) == 1
                    text = strcat(text, names.v{i}, " ");
                end
            end
            text = strcat(text, "\n");

            % Tasks (with completion times)
            text = strcat(text, "- Tasks:");
            Gam = self.task_times(z);
            if sum(z(self.dv.n_e + self.dv.n + (1:self.dv.n_t))) == 0
                text = strcat(text, " none\n");
            else
                text = strcat(text, "\n");
                for i = 1:self.dv.n_t
                    if z(self.dv.n_e + self.dv.n + i) == 1
                        [location, task] = find(self.idx_Evt == i);
                        text = strcat(text, "  > Task ", names.t{task}, " at location ", names.v{location}, ...
                            sprintf(" (completed at time %.4g)", Gam(location, task)), "\n");
                    end
                end
            end

            % Ordered list of visited locations
            text = strcat(text, "- Path:");
            path = self.get_path(z);
            for i = 1:length(path)
                text = strcat(text, " ", names.v{path(i)});
                if i < length(path); text = strcat(text, " ->"); end
            end
            text = strcat(text, "\n----------------\n");

        end

        function idx = get_end_node(self, z)
            % Returns the index of the graph Gv corresponding to the end-node of the agent for
            % the partial solution @z
            idx = find(z(self.dv.n_e + (1:self.dv.n)) == 1, 1);
            if isempty(idx); idx = 0; end
        end

        function idx = get_movement_edges(self, z)
            % Returns the edge indices of the graph Gv corresponding to the movements
            % performed by the agent as indicated in the partial solution @z
            idx = find(z(1:self.dv.n_e) == 1);
        end

        function [l_idx, t_idx] = get_task_idx(self, z)
            % Returns the indices of the location-tasks of the agent corresponding to
            % the partial solution @z
            done = find(z(self.dv.n_e + self.dv.n + (1:self.dv.n_t)) == 1);
            l_idx = zeros(1, length(done));
            t_idx = zeros(1, length(done));
            for j = 1:length(done)
                [l_idx(j), t_idx(j)] = find(self.idx_Evt == done(j));
            end
        end

    end

    methods (Static)

        function sets = short_cycles(adj, max_sets)
            % Returns the sets of 3 or 4 locations that form a cycle of the undirected graph with adjacency
            % matrix @adj (directed adjacency; the direction is ignored), at most @max_sets of them (the
            % triangles first). Each row of @sets contains the (sorted) locations of a set, padded with 0.
            n = size(adj, 1);
            U = (adj + adj') > 0;
            U(1:(n+1):end) = false;
            nb = cell(1, n);
            for v = 1:n
                nb{v} = find(U(:, v))';
            end
            tri = zeros(0, 3);
            quad = zeros(0, 4);
            for a = 1:n
                na = nb{a}(nb{a} > a); % a is the smallest location of the set
                for p = 1:length(na)
                    for q = (p+1):length(na)
                        b = na(p); d = na(q);
                        if U(b, d)
                            tri(end+1, :) = [a, b, d]; %#ok<AGROW>
                        end
                        common = nb{b}(U(nb{b}, d)); % Locations adjacent to b and d (cycle a-b-c-d-a)
                        common = common(common > a);
                        for c = common
                            quad(end+1, :) = sort([a, b, c, d]); %#ok<AGROW>
                        end
                    end
                end
                if size(tri, 1) + size(quad, 1) > 4*max_sets; break; end
            end
            quad = unique(quad, 'rows');
            sets = [tri, zeros(size(tri, 1), 1); quad];
            sets = sets(1:min(size(sets, 1), max_sets), :);
        end

    end

end
