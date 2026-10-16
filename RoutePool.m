%% A class for storing the routes of the agents (route pool) and building the route selection problem
%
% Used by primal_recovery_MRTA. For each agent, the pool contains its routes (see RouteGraph), their
% cost (8) without the task costs, and the pending tasks that the agent can perform on each route.
%
% If 'prefixes' is true, adding a route also adds its prefixes that end at its first location or at a
% location where the agent can perform a pending task, since the agent can stop at any of them. The
% route itself is only added if it ends at one of those locations (otherwise, its longest such prefix
% performs the same tasks at a lower cost).
%

classdef RoutePool < handle

    properties
        m % Number of agents
        m_z % Number of pending tasks
        graph % RouteGraph used to compute the cost of the routes
        task_v % m_z x 1 vector with the location of each pending task
        can % m_z x m logical matrix: can(j, i) is true if agent i can perform pending task j
        task_loc % m x n logical matrix: locations where each agent can perform at least one pending task
        prefixes % Add the prefixes of the routes (see above)
        paths % 1 x m cell of cell arrays with the routes of each agent
        E % 1 x m cell of column vectors with the cost of each route (without the task costs)
        cov % 1 x m cell of cell arrays with the pending tasks that the agent can perform on each route
        key_list % 1 x m cell of cell arrays with the key of each route
        hash % 1 x m cell of column vectors with a hash of each route (see route_hashes), used to find routes
    end

    methods

        function obj = RoutePool(graph, task_v, can, prefixes)
            % RoutePool(graph, task_v, can, prefixes)
            %   @graph: RouteGraph of the location graph
            %   @task_v: Location of each pending task
            %   @can: m_z x m logical matrix: can(j, i) is true if agent i can perform pending task j
            %   @prefixes: Add the prefixes of the routes (see above)
            obj.graph = graph;
            obj.task_v = reshape(task_v, [], 1);
            obj.can = logical(can);
            [obj.m_z, obj.m] = size(obj.can);
            obj.prefixes = logical(prefixes);
            obj.task_loc = false(obj.m, graph.n);
            for i = 1:obj.m
                obj.task_loc(i, obj.task_v(obj.can(:, i))) = true;
            end
            obj.paths = cell(1, obj.m);
            obj.E = cell(1, obj.m);
            obj.cov = cell(1, obj.m);
            obj.key_list = cell(1, obj.m);
            obj.hash = cell(1, obj.m);
            for i = 1:obj.m
                obj.paths{i} = cell(1, 0);
                obj.E{i} = zeros(0, 1);
                obj.cov{i} = cell(1, 0);
                obj.key_list{i} = cell(1, 0);
                obj.hash{i} = zeros(0, 1);
            end
        end

        function n_new = add(self, i, route)
            % Adds the route @route (row vector of locations) of agent @i (and its prefixes, see above).
            % Returns the number of routes that were not in the pool.
            route = reshape(route, 1, []);
            if self.prefixes
                last = [1, find(self.task_loc(i, route(2:end))) + 1];
            else
                last = length(route);
            end
            n_new = 0;
            h = RoutePool.route_hashes(route); % h(q): hash of route(1:q)
            % The prefixes of the routes of the pool are also in the pool, so the prefixes are added from
            % the longest one until one of them is already in the pool
            for q = fliplr(last)
                r = route(1:q);
                if self.find_route(i, r, h(q)) > 0
                    break;
                end
                idx = length(self.paths{i}) + 1;
                self.hash{i}(idx, 1) = h(q);
                self.key_list{i}{idx} = RoutePool.route_key(r);
                self.paths{i}{idx} = r;
                self.E{i}(idx, 1) = self.graph.route_cost(r);
                visited = false(1, self.graph.n);
                visited(r) = true;
                on_route = visited(self.task_v);
                self.cov{i}{idx} = reshape(find(self.can(:, i) & on_route(:)), [], 1);
                n_new = n_new + 1;
            end
        end

        function idx = find_route(self, i, route, h)
            % Returns the index of the route @route (row vector) of agent @i in the pool (0 if it is not in
            % the pool). @h is its hash (see route_hashes)
            idx = 0;
            for c = reshape(find(self.hash{i} == h), 1, [])
                if isequal(self.paths{i}{c}, route)
                    idx = c;
                    return
                end
            end
        end

        function n = num_routes(self)
            % Total number of routes in the pool
            n = sum(cellfun(@length, self.paths));
        end

        function [miqp, col_agent, col_route, offset] = selection_problem(self, rho)
            % Route selection problem (18) in the form used in the article when the task costs do not
            % depend on the robot (Section VI-A), with a binary variable u_c for each route c of the pool:
            %   min  sum_c E_c*u_c + sum_j rho_{t_j}
            %   s.t. sum_{c of agent i} u_c = 1,                         for each agent i
            %        sum_{c visiting v_j, agent of c can do j} u_c >= 1, for each pending task j
            % The task costs sum_j rho_{t_j} (@offset) do not depend on the assignment, which is obtained
            % afterwards with assign_tasks. The routes of each agent are ordered by their key.
            % Returns the MIQPholder @miqp, the agent and index (in self.paths) of the route of each
            % variable, and the constant term @offset of the cost.
            %   @rho: Cost rho_{t_j} of each pending task
            m = self.m;
            m_z = self.m_z;
            n_y = self.num_routes();
            col_agent = zeros(n_y, 1);
            col_route = zeros(n_y, 1);
            E_y = zeros(n_y, 1);
            cov_rows = cell(m, 1);
            cov_cols = cell(m, 1);
            cnt = 0;
            for i = 1:m
                [~, ord] = sort(self.key_list{i});
                ord = ord(:);
                idx = cnt + (1:length(ord))';
                col_agent(idx) = i;
                col_route(idx) = ord;
                E_y(idx) = self.E{i}(ord);
                cov_i = self.cov{i}(ord);
                lens = cellfun(@length, cov_i);
                cov_rows{i} = zeros(0, 1);
                cov_cols{i} = zeros(0, 1);
                if any(lens > 0)
                    cov_rows{i} = reshape(vertcat(cov_i{:}), [], 1);
                    cov_cols{i} = reshape(repelem(idx(:), lens(:)), [], 1); % repelem returns a row if idx is a scalar
                end
                cnt = cnt + length(ord);
            end
            offset = sum(rho);
            C = sparse(col_agent, (1:n_y)', 1, m, n_y);
            b = ones(m, 1);
            A = sparse(vertcat(cov_rows{:}), vertcat(cov_cols{:}), 1, m_z, n_y);
            lb = ones(m_z, 1);
            ub = inf(m_z, 1);
            miqp = MIQPholder([], E_y, A, lb, ub, C, b, [], [], [], zeros(n_y, 1), ones(n_y, 1), (1:n_y)');
        end

        function assign = assign_tasks(self, route_idx)
            % Returns the assignment of the pending tasks (assign(j) = agent that performs task j, 0 if no
            % agent can) for the routes route_idx(i) of the agents (indices in self.paths): each task goes
            % to the agent whose route arrives first at its location (ties: smallest index)
            assign = zeros(self.m_z, 1);
            best = inf(self.m_z, 1);
            for i = 1:self.m
                route = self.paths{i}{route_idx(i)};
                arrival = inf(1, self.graph.n); % Arrival time of the route at each location
                arrival(route) = [0, cumsum(self.graph.W(route(1:end-1) + (route(2:end) - 1)*self.graph.n))];
                S = self.cov{i}{route_idx(i)};
                t = reshape(arrival(self.task_v(S)), [], 1);
                upd = t < best(S);
                assign(S(upd)) = i;
                best(S(upd)) = t(upd);
            end
        end

    end

    methods (Static)

        function key = route_key(route)
            % Key of a route (same as in the route pool of heur_distLag_MRTA)
            key = sprintf('%d,', route);
        end

        function h = route_hashes(route)
            % Returns h(q) = hash of route(1:q) for each q (polynomial hash modulo the prime 2^31 - 1, exact
            % in double precision). Routes with the same hash are compared element by element.
            p = 2147483647;
            h = zeros(1, length(route));
            v = 0;
            for q = 1:length(route)
                v = mod(v*131 + route(q), p);
                h(q) = v;
            end
        end

    end

end
