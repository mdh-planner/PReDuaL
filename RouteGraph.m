%% A class with the graph tools used by the route search of the primal recovery
%
% A route of an agent is a simple path of the location graph (row vector of location indices) that
% starts at the location of the agent, i.e., the path of a solution of its MILP problem. The cost (8)
% of a route without the task costs is
%   E = w_{i,i} + sum of the weights w_{i,j} of its edges + epsilon_tau*(sum of its accumulated travel times).
%
% neighbors(), segment_removals(), tail_transfers() and segment_transfers() generate the routes (a)-(d)
% of Section VI-B from a given route. The locations are connected by shortest paths that do not visit
% the rest of the route, so the new routes are also simple paths. The shortest paths are computed once
% with Floyd-Warshall, and Dijkstra is only used when a shortest path visits a location to be avoided.
%

classdef RouteGraph < handle

    properties
        n % Number of locations
        W % n x n matrix with the weight of each edge (NaN if there is no edge)
        out_nb % 1 x n cell: successors of each location (row vectors)
        out_w % 1 x n cell: weights of the edges to the successors
        in_nb % 1 x n cell: predecessors of each location (row vectors)
        in_w % 1 x n cell: weights of the edges from the predecessors
        D % n x n matrix with the length of the shortest path between each pair of locations
        NX % n x n matrix: NX(u, v) is the next location of the shortest path from u to v (0 if none)
        cost_stop % Cost w_{i,i} of stopping at a location
        cost_times % Cost epsilon_tau of the time variables tau
        SP % n x n cell: SP{u, v} is the shortest path from u to v (filled by sp_path when needed)
    end

    methods

        function obj = RouteGraph(Gv, param)
            % RouteGraph(Gv, param)
            %   @Gv: Location graph (digraph)
            %   @param: Structure with the fields cost_stop and cost_times (see MRTA.default_param)
            n_v = numnodes(Gv);
            obj.n = n_v;
            [s, t] = findedge(Gv);
            s = reshape(s, 1, []);
            t = reshape(t, 1, []);
            w = reshape(Gv.Edges.Weight, 1, []);
            obj.W = NaN(n_v);
            obj.W(s + (t - 1)*n_v) = w;
            obj.out_nb = cell(1, n_v);
            obj.out_w = cell(1, n_v);
            obj.in_nb = cell(1, n_v);
            obj.in_w = cell(1, n_v);
            for v = 1:n_v
                obj.out_nb{v} = t(s == v);
                obj.out_w{v} = w(s == v);
                obj.in_nb{v} = s(t == v);
                obj.in_w{v} = w(t == v);
            end
            obj.cost_stop = param.cost_stop;
            obj.cost_times = param.cost_times;

            % All-pairs shortest paths (Floyd-Warshall) with next-location matrix
            dist = inf(n_v);
            dist(1:(n_v + 1):end) = 0;
            nx = zeros(n_v);
            for e = 1:length(w)
                if s(e) ~= t(e) && w(e) < dist(s(e), t(e))
                    dist(s(e), t(e)) = w(e);
                    nx(s(e), t(e)) = t(e);
                end
            end
            for k = 1:n_v
                alt = dist(:, k) + dist(k, :);
                upd = alt < dist;
                if any(upd(:))
                    dist(upd) = alt(upd);
                    nx_k = repmat(nx(:, k), 1, n_v);
                    nx(upd) = nx_k(upd);
                end
            end
            obj.D = dist;
            obj.NX = nx;
            obj.SP = cell(n_v);
        end

        function E = route_cost(self, route)
            % Cost (8) of the route @route without the task costs
            E = self.cost_stop;
            if length(route) > 1
                w = self.W(route(1:end-1) + (route(2:end) - 1)*self.n);
                E = E + sum(w) + self.cost_times*sum(cumsum(w));
            end
        end

        function p = sp_path(self, u, v)
            % Shortest path [u, ..., v] from location u to location v (empty if there is none)
            if u == v
                p = u;
                return
            end
            if ~isfinite(self.D(u, v))
                p = [];
                return
            end
            p = self.SP{u, v};
            if ~isempty(p); return; end
            p = zeros(1, self.n);
            p(1) = u;
            len = 1;
            while p(len) ~= v
                len = len + 1;
                if len > self.n
                    error("RouteGraph:sp_path", "Invalid shortest path from %d to %d", u, v);
                end
                p(len) = self.NX(p(len-1), v);
            end
            p = p(1:len);
            self.SP{u, v} = p;
        end

        function p = path_avoiding(self, u, v, forbidden)
            % Shortest path [u, ..., v] from location u to location v that does not visit the locations
            % with forbidden(.) = true (except u). Empty if there is none.
            p = self.sp_path(u, v);
            if isempty(p) || ~any(forbidden(p(2:end)))
                return
            end
            [dist, pred] = self.shortest_paths(u, forbidden, false, v);
            if isfinite(dist(v))
                p = RouteGraph.trace(pred, u, v);
            else
                p = [];
            end
        end

        function [dist, pred] = shortest_paths(self, src, forbidden, reverse, targets)
            % [dist, pred] = shortest_paths(src, forbidden, reverse, [targets])
            % Dijkstra's algorithm in the graph without the locations v with forbidden(v) = true (the
            % location @src is never forbidden). If @reverse is false, dist(v) is the length of the
            % shortest path from @src to v and pred(v) the previous location of v in that path. If
            % @reverse is true, dist(v) is the length of the shortest path from v to @src and pred(v)
            % the next location of v in that path. dist(v) = Inf and pred(v) = 0 if there is no path.
            % If @targets (vector of locations) is given, the algorithm stops when the shortest paths of all
            % the targets are known: dist and pred are then only valid for the targets and the locations
            % closer to @src than the farthest target (the paths to the targets are the same as without
            % @targets).
            n_v = self.n;
            dist = inf(1, n_v);
            pred = zeros(1, n_v);
            closed = reshape(logical(forbidden), 1, n_v);
            closed(src) = false;
            dist(src) = 0;
            is_target = false(1, n_v);
            if nargin >= 5
                is_target(targets) = true;
                is_target(closed) = false; % Forbidden locations are never reached
            else
                is_target(:) = true;
            end
            num_left = nnz(is_target);
            if reverse
                nbs = self.in_nb;
                ws = self.in_w;
            else
                nbs = self.out_nb;
                ws = self.out_w;
            end
            for it = 1:n_v
                d = dist;
                d(closed) = inf;
                [d_u, u] = min(d);
                if ~isfinite(d_u); break; end
                if is_target(u)
                    num_left = num_left - 1;
                    if num_left == 0; break; end
                end
                closed(u) = true;
                v = nbs{u};
                d_v = d_u + ws{u};
                upd = ~closed(v) & d_v < dist(v);
                dist(v(upd)) = d_v(upd);
                pred(v(upd)) = u;
            end
        end

        function routes = neighbors(self, route, anchor, cand)
            % routes = neighbors(route, anchor, cand)
            % Returns a cell array with the following routes obtained from @route (routes (a) and (b)):
            %   - Removal: for each anchor location of the route (except its first location), the route
            %     without it. Its predecessor and successor anchor locations are connected by a shortest
            %     path that does not visit the rest of the route (if it is the last anchor location, the
            %     route ends at the previous one).
            %   - Insertion: for each location v with cand(v) = true that is not in the route, the
            %     cheapest route obtained by visiting v between two consecutive anchor locations, or
            %     after the last one, using shortest paths that do not visit the rest of the route.
            %   @route: Row vector with the route (simple path)
            %   @anchor: 1 x n logical vector with the locations that must be kept (those where the agent
            %            performs tasks). The first location of the route is always kept.
            %   @cand: 1 x n logical vector with the locations that can be inserted
            route = reshape(route, 1, []);
            n_v = self.n;
            pos = [1, find(anchor(route(2:end))) + 1]; % Positions of the kept locations
            routes = cell(1, 0);

            % Removal
            for g = 2:length(pos)
                if g == length(pos)
                    routes{end+1} = route(1:pos(g-1)); %#ok<AGROW>
                else
                    pre = route(1:pos(g-1));
                    suf = route(pos(g+1):end);
                    forbidden = false(1, n_v);
                    forbidden([pre(1:end-1), suf(2:end)]) = true;
                    leg = self.path_avoiding(pre(end), suf(1), forbidden);
                    if ~isempty(leg)
                        routes{end+1} = [pre, leg(2:end-1), suf]; %#ok<AGROW>
                    end
                end
            end

            % Insertion
            on_route = false(1, n_v);
            on_route(route) = true;
            vs = find(reshape(cand, 1, n_v) & ~on_route);
            if isempty(vs); return; end
            best_E = inf(1, length(vs));
            best_route = cell(1, length(vs));
            for g = 1:length(pos)
                a = route(pos(g));
                pre = route(1:pos(g));
                forbidden = false(1, n_v);
                forbidden(pre(1:end-1)) = true;
                if g == length(pos)
                    % After the last anchor location
                    dist = []; % Shortest paths from a avoiding the route (computed only if needed)
                    for q = 1:length(vs)
                        leg = self.sp_path(a, vs(q));
                        if isempty(leg); continue; end
                        if any(forbidden(leg))
                            if isempty(dist)
                                [dist, pred] = self.shortest_paths(a, forbidden, false, vs);
                            end
                            if ~isfinite(dist(vs(q))); continue; end
                            leg = RouteGraph.trace(pred, a, vs(q));
                        end
                        r = [pre, leg(2:end)];
                        E = self.route_cost(r);
                        if E < best_E(q)
                            best_E(q) = E;
                            best_route{q} = r;
                        end
                    end
                else
                    % Between the anchor locations a and b
                    suf = route(pos(g+1):end);
                    b = suf(1);
                    forbidden(suf(2:end)) = true;
                    f_a = forbidden; % Locations that the path from a to v cannot visit
                    f_a(b) = true;
                    f_b = forbidden; % Locations that the path from v to b cannot visit
                    f_b(a) = true;
                    d_1 = [];
                    d_2 = [];
                    for q = 1:length(vs)
                        v = vs(q);
                        leg_1 = self.sp_path(a, v); % a -> v
                        leg_2 = self.sp_path(v, b); % v -> b
                        if isempty(leg_1) || isempty(leg_2); continue; end
                        if any(f_a(leg_1))
                            if isempty(d_1)
                                [d_1, p_1] = self.shortest_paths(a, f_a, false, vs);
                            end
                            if ~isfinite(d_1(v)); continue; end
                            leg_1 = RouteGraph.trace(p_1, a, v);
                        end
                        if any(f_b(leg_2))
                            if isempty(d_2)
                                [d_2, p_2] = self.shortest_paths(b, f_b, true, vs);
                            end
                            if ~isfinite(d_2(v)); continue; end
                            leg_2 = RouteGraph.trace_reverse(p_2, b, v);
                        end
                        if any(ismember(leg_1(1:end-1), leg_2))
                            % The paths intersect: path from v to b that does not visit leg_1
                            f_v = forbidden;
                            f_v(leg_1(1:end-1)) = true;
                            leg_2 = self.path_avoiding(v, b, f_v);
                            if isempty(leg_2); continue; end
                        end
                        r = [pre, leg_1(2:end), leg_2(2:end-1), suf];
                        E = self.route_cost(r);
                        if E < best_E(q)
                            best_E(q) = E;
                            best_route{q} = r;
                        end
                    end
                end
            end
            routes = [routes, best_route(isfinite(best_E))];
        end

        function routes = segment_removals(self, route, anchor, max_len)
            % routes = segment_removals(route, anchor, max_len)
            % Returns a cell array with the routes obtained from @route by removing 2 to @max_len
            % consecutive anchor locations (see neighbors), other than the last one. Their predecessor and
            % successor anchor locations are connected by a shortest path that does not visit the rest of
            % the route. Removing the last anchor locations gives a prefix of the route, which is already
            % in the pool (see RoutePool).
            route = reshape(route, 1, []);
            n_v = self.n;
            pos = [1, find(anchor(route(2:end))) + 1];
            n_pos = length(pos);
            routes = cell(1, 0);
            for g_1 = 2:n_pos
                for g_2 = (g_1 + 1):min(g_1 + max_len - 1, n_pos - 1)
                    pre = route(1:pos(g_1-1));
                    suf = route(pos(g_2+1):end);
                    forbidden = false(1, n_v);
                    forbidden([pre(1:end-1), suf(2:end)]) = true;
                    leg = self.path_avoiding(pre(end), suf(1), forbidden);
                    if ~isempty(leg)
                        routes{end+1} = [pre, leg(2:end-1), suf]; %#ok<AGROW>
                    end
                end
            end
        end

        function routes = tail_transfers(self, route, anchor, others, can_do)
            % routes = tail_transfers(route, anchor, others, can_do)
            % Routes (c): returns a cell array with the routes in which the agent travels the beginning of its route
            % @route (up to its first location or up to one of its anchor locations) followed by the end of
            % the route of another agent (from one of the anchor locations of that route), connected by a
            % shortest path that does not visit the rest of the new route. Combined with the routes in
            % which the other agent stops earlier (prefixes, see RoutePool), they allow the route
            % selection to transfer several tasks at once from one agent to another, or to exchange the
            % ends of two routes.
            %   @route, @anchor: Route and anchor locations of the agent (see neighbors)
            %   @others: Cell array where others{k} = {route_k, anchor_k} for each other agent
            %   @can_do: 1 x n logical vector with the locations where the agent can perform tasks.
            %            Only the ends of routes that visit at least one of these locations are used.
            route = reshape(route, 1, []);
            pos = [1, find(anchor(route(2:end))) + 1];
            pieces = cell(1, 0);
            for k = 1:length(others)
                r_k = reshape(others{k}{1}, 1, []);
                anchor_k = others{k}{2};
                for g = find(anchor_k(r_k(2:end))) + 1
                    tail = r_k(g:end);
                    if any(can_do(tail))
                        pieces{end+1} = tail; %#ok<AGROW>
                    end
                end
            end
            routes = self.append_pieces(route, pos, pieces);
        end

        function routes = segment_transfers(self, route, anchor, others, can_do, max_len)
            % routes = segment_transfers(route, anchor, others, can_do, max_len)
            % Routes (d): returns a cell array with the routes in which the agent travels its route @route (up to its
            % last anchor location) followed by a part of the route of another agent, from one of its anchor
            % locations to another one (1 to @max_len anchor locations, without the last one: see
            % tail_transfers), connected by a shortest path that does not visit the rest of the new route.
            % Combined with the routes of the other agent without that part (see segment_removals), they
            % allow the route selection to transfer groups of consecutive tasks between agents.
            %   @route, @anchor, @others, @can_do: See tail_transfers
            route = reshape(route, 1, []);
            pos = [1, find(anchor(route(2:end))) + 1];
            pieces = cell(1, 0);
            for k = 1:length(others)
                r_k = reshape(others{k}{1}, 1, []);
                anchor_k = others{k}{2};
                p_k = find(anchor_k(r_k(2:end))) + 1;
                n_k = length(p_k);
                for x = 1:n_k
                    for y = x:min(x + max_len - 1, n_k - 1)
                        seg = r_k(p_k(x):p_k(y));
                        if any(can_do(seg))
                            pieces{end+1} = seg; %#ok<AGROW>
                        end
                    end
                end
            end
            routes = self.append_pieces(route, pos(end), pieces);
        end

    end

    methods (Access = private)

        function routes = append_pieces(self, route, cuts, pieces)
            % Returns a cell array with the routes route(1:h) + (shortest path) + piece, for each position h
            % in @cuts and each path piece{q} (row vector) that does not visit route(1:h), except for its
            % first location, which can be route(h) (the route is then route(1:h) + piece(2:end), e.g., if
            % two tasks at the same location are performed by different agents). The shortest path from
            % route(h) to piece(1) does not visit route(1:h-1) and piece(2:end).
            n_v = self.n;
            routes = cell(1, 0);
            if isempty(pieces); return; end
            for h = reshape(cuts, 1, [])
                pre = route(1:h);
                a = pre(end);
                in_pre = false(1, n_v);
                in_pre(pre) = true;
                forbidden = in_pre;
                forbidden(a) = false;
                dist = []; % Shortest paths from a avoiding the beginning of the route (computed only if needed)
                firsts = cellfun(@(piece) piece(1), pieces); % First location of each piece
                for q = 1:length(pieces)
                    piece = pieces{q};
                    if any(in_pre(piece(2:end))); continue; end % The piece visits the beginning of the route
                    if piece(1) == a
                        routes{end+1} = [pre, piece(2:end)]; %#ok<AGROW>
                        continue
                    end
                    if in_pre(piece(1)); continue; end
                    f_p = forbidden;
                    f_p(piece(2:end)) = true;
                    leg = self.sp_path(a, piece(1));
                    if isempty(leg); continue; end
                    if any(f_p(leg))
                        if isempty(dist)
                            [dist, pred] = self.shortest_paths(a, forbidden, false, firsts);
                        end
                        if ~isfinite(dist(piece(1))); continue; end
                        leg = RouteGraph.trace(pred, a, piece(1));
                        if any(f_p(leg))
                            leg = self.path_avoiding(a, piece(1), f_p);
                            if isempty(leg); continue; end
                        end
                    end
                    routes{end+1} = [pre, leg(2:end), piece(2:end)]; %#ok<AGROW>
                end
            end
        end

    end

    methods (Static)

        function p = trace(pred, src, dst)
            % Returns the path [src, ..., dst] of a shortest-path tree from src (pred computed by
            % shortest_paths with reverse = false)
            p = zeros(1, length(pred));
            p(1) = dst;
            len = 1;
            while p(len) ~= src
                len = len + 1;
                if len > length(pred) || pred(p(len-1)) == 0
                    error("RouteGraph:trace", "Location %d cannot be reached", dst);
                end
                p(len) = pred(p(len-1));
            end
            p = p(len:-1:1);
        end

        function p = trace_reverse(pred, dst, src)
            % Returns the path [src, ..., dst] of a shortest-path tree to dst (pred computed by
            % shortest_paths with reverse = true)
            p = zeros(1, length(pred));
            p(1) = src;
            len = 1;
            while p(len) ~= dst
                len = len + 1;
                if len > length(pred) || pred(p(len-1)) == 0
                    error("RouteGraph:trace", "Location %d cannot reach location %d", src, dst);
                end
                p(len) = pred(p(len-1));
            end
            p = p(1:len);
        end

    end

end
