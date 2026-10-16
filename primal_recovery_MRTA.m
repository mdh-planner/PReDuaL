%% primal_recovery_MRTA - Primal recovery of PReDuaL (Section VI and Algorithm 3 of the article)
%
% Builds a feasible solution of the central MILP (9) from the routes computed by the agents during
% Algorithm 1 (route pool). Each round of the recovery does the following:
%   1. Route selection: solve problem (18), which picks one route of the pool per agent and assigns
%      each pending task to an agent whose route visits its location (see RoutePool.selection_problem).
%   2. Route search (opt.route_search): generate new routes from the selected ones without solving any
%      MILP problem (routes (a)-(d) of Section VI-B, see RouteGraph), add them to the pool and solve (18)
%      again, while its cost decreases.
%   3. Re-routing: each agent computes its optimal route for the tasks assigned to it.
%   4. Local search (opt.local_search): move single tasks between agents while the cost decreases.
%   5. Pairwise re-optimization (opt.pair_reopt): redistribute the tasks of each pair of agents by
%      solving (9) restricted to the pair.
% The routes computed in steps 3-5 are added to the pool, and the rounds are repeated while the cost
% decreases (at most opt.rounds times). If (18) is infeasible, the last iterate of Algorithm 1 is
% repaired instead: a task done by several agents is kept by the one that completes it first, and a
% task done by none is inserted where it is cheapest.
%
% Optionally (opt.column_generation), the rounds are followed by column generation on the LP
% relaxation of (18): the agents price new routes with their MILP problems using the dual variables of
% the LP, which also gives a Lagrangian lower bound of the optimal cost of (9). If new routes are
% found, the rounds are repeated with the enlarged pool.
%
% INPUTS:
%   @mrta: Instance of MRTA with the agents built by build_distMILP
%   @idx_A: Structure with the location (.v) and task-type (.t) of each pending task (rows of (3))
%   @pool: 1 x m cell of containers.Map (key -> path) with the routes computed by each agent
%   @Gam: m_z x m matrix with the completion times Gamma_i(t_j) of the last iterate (Inf if not done)
%   @opt: Structure with fields
%       .rounds: Maximum number of rounds
%       .route_search: Apply the route search (default false)
%       .route_dp: Compute routes with RouteEvaluator.route_dp when possible (default true)
%       .rs_max_iter: Maximum number of repetitions of the route search in each round (default 50)
%       .rs_max_segment: Maximum number of consecutive task locations of routes (a) and (d) (default 3)
%       .local_search: Apply the local search
%       .ls_max_pass: Maximum number of passes of the local search
%       .pair_reopt: Apply the pairwise re-optimization
%       .column_generation: Apply the column generation
%       .cg_max_iter: Maximum number of column generation rounds
%       .cg_smoothing: Dual smoothing parameter in [0, 1) (0 for no smoothing)
%       .cg_tol: Relative tolerance used to stop the column generation
%       .solver_args: Cell array with name/value options for MIQPholder.solve
%       .agent_solver: (optional) AgentSolver for the MILP problems of the agents (sequential if not given)
%
% OUTPUTS:
%   @z: Feasible solution of (9) (empty if none was found)
%   @V: Cost (8) of @z (Inf if no feasible solution was found)
%   @rec: Structure with information:
%       .LB: Lower bound of the optimal cost of (9) (-Inf if not computed)
%       .gap: Relative optimality gap (V - LB)/LB
%       .V_phase1: Cost before the column generation
%       .V_select, .V_search, .V_reroute, .V_ls, .V_pair: Cost after each step of each round
%       .used_repair: true if the last iterate of Algorithm 1 was repaired
%       .cg: Information of the column generation
%       .num_solves: Number of MILP problems solved by the agents
%       .num_dp: Number of routes computed by RouteEvaluator.route_dp
%       .num_selections: Number of route selection problems solved
%       .pool_size: Number of routes in the pool at the end
%       .run_time_search: Time spent generating the routes of the route search
%       .run_time_par_total: Computation time, as in Section VI-C
%

function [z, V, rec] = primal_recovery_MRTA(mrta, idx_A, pool, Gam, opt)

    start_time = tic;
    m = mrta.m;
    agent_solver = [];
    if isfield(opt, 'agent_solver'); agent_solver = opt.agent_solver; end
    if ~isfield(opt, 'route_search'); opt.route_search = false; end
    if ~isfield(opt, 'rs_max_iter'); opt.rs_max_iter = 50; end
    if ~isfield(opt, 'rs_max_segment'); opt.rs_max_segment = 3; end
    ev = RouteEvaluator(mrta, idx_A, opt.solver_args, agent_solver);

    rec.V_select = zeros(1, 0);
    rec.V_search = zeros(1, 0);
    rec.V_reroute = zeros(1, 0);
    rec.V_ls = zeros(1, 0);
    rec.V_pair = zeros(1, 0);
    rec.used_repair = false;
    rec.run_time_coord = 0; % Solver time of the route selection problems (done by the coordinator)
    rec.wall_time_coord = 0;
    rec.run_time_search = 0;
    rec.num_selections = 0;

    %% Route pool: routes computed by the agents in Algorithm 1
    graph = RouteGraph(mrta.Gv, mrta.agent{1}.param);
    ev.graph = graph;
    if isfield(opt, 'route_dp'); ev.use_dp = opt.route_dp; end
    rp = RoutePool(graph, idx_A.v, ev.idx_loc > 0, opt.route_search);
    for i = 1:m
        paths = values(pool{i});
        for p = 1:length(paths)
            rp.add(i, paths{p});
        end
        if opt.route_search
            % Shortest path from the agent to each location where it can perform a pending task
            t_start = tic;
            seeds = graph.neighbors(mrta.agent{i}.loc, false(1, graph.n), rp.task_loc(i, :));
            for p = 1:length(seeds)
                rp.add(i, seeds{p});
            end
            rec.run_time_search = rec.run_time_search + toc(t_start);
        end
    end

    %% Rounds of route selection, route search, re-routing and local search
    searched = containers.Map('KeyType', 'char', 'ValueType', 'logical'); % See route_search
    [assign, V, rec] = improve(mrta, rp, graph, Gam, ev, opt, [], inf, rec, searched);
    rec.V_phase1 = V;

    %% Column generation (new routes and a lower bound)
    rec.LB = -inf;
    rec.cg = struct('iters', 0, 'num_new', 0, 'V_lp', NaN, 'converged', false);
    if opt.column_generation && ~isempty(assign)
        [LB, cg, rec] = column_generation(rp, ev, V, opt, rec);
        rec.LB = LB;
        rec.cg = cg;
        if rec.cg.num_new > 0 && V - rec.LB > opt.cg_tol*max(1, abs(V))
            [assign, V, rec] = improve(mrta, rp, graph, [], ev, opt, assign, V, rec, searched);
        end
    end
    rec.V_phase2 = V;

    %% Construct the solution
    if isempty(assign)
        z = [];
        V = inf;
    else
        z = zeros(mrta.nz, 1);
        V = 0;
        for i = 1:m
            ag = mrta.agent{i};
            [V_i, z_i] = ev.cost(i, find(assign == i));
            z(ag.dv.idx + (1:ag.dv.nz)) = z_i;
            V = V + V_i;
        end
    end

    %% Information
    rec.assign = assign;
    rec.gap = NaN;
    if isfinite(rec.LB) && isfinite(V)
        rec.gap = (V - rec.LB)/max(abs(rec.LB), 1e-12);
    end
    rec.num_solves = ev.num_solves;
    rec.num_dp = ev.num_dp;
    rec.pool_size = rp.num_routes();
    rec.run_time_all = toc(start_time);
    rec.run_time_par = rec.run_time_coord + ev.run_time_par;
    rec.run_time_seq = rec.run_time_coord + ev.run_time_seq;
    rec.run_time_overhead = max(rec.run_time_all - rec.wall_time_coord - ev.run_time_wall, 0);
    rec.run_time_par_total = rec.run_time_par + rec.run_time_overhead;

end

function [assign_best, V_best, rec] = improve(mrta, rp, graph, Gam, ev, opt, assign_best, V_best, rec, searched)
    % Rounds of the primal recovery, repeated while the cost decreases. Starts from the incumbent
    % assignment @assign_best with cost @V_best (empty and Inf if there is none).
    m = mrta.m;
    for r_idx = 1:opt.rounds

        % Route selection
        [sel, ok, rec] = select_routes(rp, ev, rec);
        if ok
            rec.V_select(end+1) = sel.V;
            if opt.route_search
                % Route search
                [sel, rec] = route_search(rp, graph, sel, ev, opt, rec, searched);
                rec.V_search(end+1) = sel.V;
            end
            assign = sel.assign;
        elseif r_idx == 1 && isempty(assign_best) && ~isempty(Gam)
            % (18) is infeasible: repair the last iterate of Algorithm 1
            assign = repair_assignment(Gam, ev);
            rec.used_repair = true;
            if isempty(assign); break; end
        else
            break;
        end

        % Re-routing: each agent computes its optimal route for its tasks. The selected route (up to the
        % location of its last task) is used to warm-start the solver.
        sets = arrayfun(@(i) find(assign == i), 1:m, 'UniformOutput', false);
        ws = cell(1, m);
        if ok
            for i = 1:m
                ws{i} = selected_route_z(mrta.agent{i}, rp.paths{i}{sel.route(i)}, sets{i}, ev, i);
            end
        end
        ev.evaluate(1:m, sets, ws);
        if ok
            % Keep the selected route if it is cheaper than the MILP solution (the solver stops at its
            % relative gap tolerance)
            for i = 1:m
                route_i = rp.paths{i}{sel.route(i)};
                V_i = graph.route_cost(route_i) + sum(ev.task_rho(sets{i}));
                if V_i < ev.cost(i, sets{i}) - 1e-9*max(1, abs(V_i))
                    ag = mrta.agent{i};
                    z_i = ag.route_to_z(route_i, ev.idx_loc(sets{i}, i));
                    if ag.miqp.is_feasible(z_i)
                        ev.store(i, sets{i}, z_i);
                    end
                end
            end
        end
        V_round = 0;
        for i = 1:m
            V_round = V_round + ev.cost(i, sets{i});
        end
        ev.end_round();
        rec.V_reroute(end+1) = V_round;

        % Local search
        if opt.local_search
            [assign, V_round] = local_search(assign, ev, opt.ls_max_pass);
            rec.V_ls(end+1) = V_round;
        end

        % Pairwise re-optimization
        if opt.pair_reopt
            [assign, V_round] = pair_reoptimization(assign, ev, opt.ls_max_pass);
            rec.V_pair(end+1) = V_round;
        end

        % Add the new routes to the pool
        num_new = add_routes_to_pool(rp, ev);

        if V_round < V_best - 1e-9
            V_best = V_round;
            assign_best = assign;
        else
            break;
        end
        if num_new == 0
            % Same pool as in the last route selection, so the next round would give the same result
            break;
        end

    end
end

function [sel, rec] = route_search(rp, graph, sel, ev, opt, rec, searched)
    % Route search: generates new routes from the selected ones (@sel, see select_routes), adds them to
    % the pool and solves the route selection problem again, while its cost decreases.
    % The new routes of agent i only depend on its route and task locations (and on those of the other
    % agent for routes (c) and (d)), so @searched keeps track of the combinations already generated.
    m = rp.m;
    for it = 1:opt.rs_max_iter
        t_start = tic;
        anchor = false(m, graph.n); % Locations where each agent performs the tasks assigned to it
        state = cell(1, m); % Key of the selected route and task locations of each agent
        for i = 1:m
            anchor(i, rp.task_v(sel.assign == i)) = true;
            state{i} = [RoutePool.route_key(rp.paths{i}{sel.route(i)}), '|', sprintf('%d,', find(anchor(i, :)))];
        end
        num_new = 0;
        for i = 1:m
            route_i = rp.paths{i}{sel.route(i)};
            routes = cell(1, 0);
            key_i = sprintf('%d#%s', i, state{i});
            if ~isKey(searched, key_i)
                routes = [graph.neighbors(route_i, anchor(i, :), rp.task_loc(i, :)), ...
                          graph.segment_removals(route_i, anchor(i, :), opt.rs_max_segment)];
                searched(key_i) = true;
            end
            others = cell(1, 0); % Other agents whose routes have not been combined with this one yet
            for k = [1:(i-1), (i+1):m]
                key_ik = sprintf('%s#%d#%s', key_i, k, state{k});
                if ~isKey(searched, key_ik)
                    others{end+1} = {rp.paths{k}{sel.route(k)}, anchor(k, :)}; %#ok<AGROW>
                    searched(key_ik) = true;
                end
            end
            if ~isempty(others)
                routes = [routes, graph.tail_transfers(route_i, anchor(i, :), others, rp.task_loc(i, :)), ...
                          graph.segment_transfers(route_i, anchor(i, :), others, rp.task_loc(i, :), opt.rs_max_segment)]; %#ok<AGROW>
            end
            for q = 1:length(routes)
                num_new = num_new + rp.add(i, routes{q});
            end
        end
        rec.run_time_search = rec.run_time_search + toc(t_start);
        if num_new == 0; break; end

        [sel_new, ok, rec] = select_routes(rp, ev, rec, sel);
        if ~ok || sel_new.V >= sel.V - 1e-9*max(1, abs(sel.V))
            break;
        end
        sel = sel_new;
    end
end

function [LB, cg, rec] = column_generation(rp, ev, UB, opt, rec)
    % Column generation on the LP relaxation of the route selection problem
    LB = -inf;
    duals_best = [];
    cg = struct('iters', 0, 'num_new', 0, 'V_lp', NaN, 'converged', false, 'hist_LB', zeros(1, 0));
    tol_lp = opt.cg_tol*max(1, abs(UB));

    for it = 1:opt.cg_max_iter
        cg.iters = it;

        % LP relaxation of the route selection problem and its dual variables
        [duals, V_lp, ok, info_lp] = selection_lp(rp, ev);
        rec.run_time_coord = rec.run_time_coord + info_lp.run_time;
        rec.wall_time_coord = rec.wall_time_coord + info_lp.wall_time;
        if ~ok; break; end
        cg.V_lp = V_lp;

        % Pricing problems of the agents: first with the smoothed duals and, if they give no new route,
        % with the duals of the LP
        candidates = {duals};
        if ~isempty(duals_best) && opt.cg_smoothing > 0
            candidates = {opt.cg_smoothing*duals_best + (1 - opt.cg_smoothing)*duals, duals};
        end
        num_new = 0;
        for q = 1:length(candidates)
            [L, valid, num_new] = price_all(candidates{q}, rp, ev);
            if valid && L > LB
                LB = L;
                duals_best = candidates{q};
            end
            if num_new > 0; break; end
        end
        cg.num_new = cg.num_new + num_new;
        cg.hist_LB(end+1) = LB;

        % Stopping criteria
        if num_new == 0 || V_lp - LB <= tol_lp
            cg.converged = true; % LP relaxation solved: its value is the best Lagrangian bound
            break;
        end
        if UB - LB <= opt.cg_tol*max(1, abs(UB))
            break; % The incumbent solution is certified
        end
    end
end

function [L, valid, num_new] = price_all(duals, rp, ev)
    % Solves the pricing problem of each agent for the dual variables @duals, adds the new routes to the
    % pool and returns the Lagrangian lower bound @L (valid if all the problems were solved to optimality)
    [~, bound, z, info] = ev.price_agents(duals);
    L = sum(duals) + sum(bound);
    valid = true;
    num_new = 0;
    for i = 1:rp.m
        valid = valid && info{i}.e_flag == 1;
        if info{i}.feasible
            route = ev.mrta.agent{i}.get_path(z{i});
            num_new = num_new + rp.add(i, route);
        end
    end
    ev.end_round();
end

function num_new = add_routes_to_pool(rp, ev)
    % Adds the routes computed by the RouteEvaluator to the pool and returns the number of new ones
    num_new = 0;
    for i = 1:rp.m
        paths = ev.get_paths(i);
        for p = 1:length(paths)
            num_new = num_new + rp.add(i, paths{p});
        end
    end
end

function [sel, ok, rec] = select_routes(rp, ev, rec, sel_prev)
    % Solves the route selection problem (18) (see RoutePool.selection_problem). Returns the structure
    % @sel with the selected route of each agent (.route, index in the pool), the assignment of the tasks
    % (.assign, see RoutePool.assign_tasks) and the optimal value (.V).
    % @sel_prev: (optional) previous solution, used to warm-start the solver (routes are only added to
    %            the pool, so it is still feasible)
    m = rp.m;
    [miqp, col_agent, col_route, offset] = rp.selection_problem(ev.task_rho);
    args = ev.solver_args;
    if nargin >= 4 && ~isempty(sel_prev) && ~isempty(sel_prev.route)
        x0 = zeros(length(col_agent), 1);
        for i = 1:m
            x0(col_agent == i & col_route == sel_prev.route(i)) = 1;
        end
        args = [args, {'warmstart', x0}];
    end
    [y, V, info] = miqp.solve(args{:});
    rec.run_time_coord = rec.run_time_coord + info.run_time;
    rec.wall_time_coord = rec.wall_time_coord + info.wall_time;
    rec.num_selections = rec.num_selections + 1;

    sel = struct('assign', [], 'route', [], 'V', V + offset);
    ok = info.feasible;
    if ok
        y = y > 0.5;
        sel.route = zeros(1, m);
        for i = 1:m
            c_i = find(y & col_agent == i, 1);
            if ~isempty(c_i)
                sel.route(i) = col_route(c_i);
            end
        end
        ok = all(sel.route > 0);
        if ok
            sel.assign = rp.assign_tasks(sel.route);
            ok = all(sel.assign > 0);
        end
    end
end

function [duals, V, ok, info] = selection_lp(rp, ev)
    % Solves the LP relaxation of the route selection problem and returns dual variables pi_j of the
    % constraints "task j is performed by exactly one agent" (one per pending task): pi_j = rho_{t_j} +
    % sigma_j, where sigma_j >= 0 is the dual variable of the covering constraint of task j (see
    % RoutePool.selection_problem). The pricing problems then use the task costs rho_{t_j} - pi_j = -sigma_j.
    [miqp, ~, ~, offset] = rp.selection_problem(ev.task_rho);
    [~, V, info] = miqp.solve(ev.solver_args{:}, 'relax', true);
    V = V + offset;
    ok = info.feasible && isfield(info, 'y_in');
    duals = [];
    if ok
        duals = ev.task_rho(:) + info.y_in(:);
    end
end

function assign = repair_assignment(Gam, ev)
    % Repairs the assignment of the last iterate of Algorithm 1 (empty if it fails)
    [m_z, m] = size(Gam);
    assign = zeros(m_z, 1);
    for j = 1:m_z
        [G_min, i] = min(Gam(j, :)); % Keep the task in the agent that completes it first
        if isfinite(G_min); assign(j) = i; end
    end
    for j = find(assign == 0)'
        % Insert the task where it is cheapest
        capable = find(ev.idx_loc(j, :) > 0);
        ev.evaluate([capable, capable], [arrayfun(@(i) find(assign == i), capable, 'UniformOutput', false), ...
                                         arrayfun(@(i) [find(assign == i); j], capable, 'UniformOutput', false)]);
        best = inf;
        i_best = 0;
        for i = 1:m
            if ~ev.can_do(j, i); continue; end
            S = find(assign == i);
            d = ev.cost(i, [S; j]) - ev.cost(i, S);
            if d < best
                best = d;
                i_best = i;
            end
        end
        ev.end_round();
        if i_best == 0
            assign = [];
            return;
        end
        assign(j) = i_best;
    end
end

function [assign, V] = local_search(assign, ev, max_pass)
    % Moves each task to the agent that reduces the cost (8) the most, until no move reduces the cost.
    % Moves that cannot reduce the cost according to the Held-Karp lower bounds are skipped.
    m_z = length(assign);
    m = size(ev.idx_loc, 2);
    for pass = 1:max_pass
        improved = false;
        for j = 1:m_z
            a = assign(j);
            S_a = find(assign == a);
            S_a_rem = S_a(S_a ~= j);
            V_a = ev.cost(a, S_a);
            delta_a_lb = ev.lower_bound(a, S_a_rem) - V_a; % Lower bound of the cost change of agent a

            % Agents for which the move can reduce the cost according to the lower bounds
            cand = zeros(1, 0);
            for b = 1:m
                if b == a || ~ev.can_do(j, b); continue; end
                S_b = find(assign == b);
                V_b = ev.cost(b, S_b);
                ins_lb = max(ev.lower_bound(b, [S_b; j]), V_b + ev.task_rho(j)); % Adding a task never reduces the route cost
                if delta_a_lb + ins_lb - V_b < -1e-9
                    cand(end+1) = b; %#ok<AGROW>
                end
            end

            % Evaluate the candidate moves
            best = -1e-9;
            b_best = 0;
            if ~isempty(cand)
                ev.evaluate([a, cand], [{S_a_rem}, arrayfun(@(b) [find(assign == b); j], cand, 'UniformOutput', false)]);
                delta_a = ev.cost(a, S_a_rem) - V_a; % Cost change of removing task j from agent a
                for b = cand
                    S_b = find(assign == b);
                    d = delta_a + ev.cost(b, [S_b; j]) - ev.cost(b, S_b);
                    if d < best
                        best = d;
                        b_best = b;
                    end
                end
            end
            ev.end_round();
            if b_best > 0
                assign(j) = b_best;
                improved = true;
            end
        end
        if ~improved; break; end
    end
    V = 0;
    for i = 1:m
        V = V + ev.cost(i, find(assign == i));
    end
end

function [assign, V] = pair_reoptimization(assign, ev, max_pass)
    % For each pair of agents (a, b), redistributes their tasks by solving problem (9) restricted to
    % agents a and b and their tasks. Repeated until no pair improves the cost.
    m = size(ev.idx_loc, 2);
    for pass = 1:max_pass
        improved = false;
        for a = 1:(m-1)
            for b = (a+1):m
                S = find(assign == a | assign == b);
                if ~any(ev.idx_loc(S, a) > 0 & ev.idx_loc(S, b) > 0)
                    continue; % No task can be exchanged between the agents
                end
                S_a = find(assign == a);
                S_b = find(assign == b);
                V_cur = ev.cost(a, S_a) + ev.cost(b, S_b);
                [ok, V_new, to_a, z_a, z_b] = solve_pair(ev, a, b, S);
                ev.end_round();
                if ok && V_new < V_cur - 1e-6*max(1, abs(V_cur))
                    assign(S(to_a)) = a;
                    assign(S(~to_a)) = b;
                    ev.store(a, S(to_a), z_a);
                    ev.store(b, S(~to_a), z_b);
                    improved = true;
                end
            end
        end
        if ~improved; break; end
    end
    V = 0;
    for i = 1:m
        V = V + ev.cost(i, find(assign == i));
    end
end

function [ok, V, to_a, z_a, z_b] = solve_pair(ev, a, b, S)
    % Solves problem (9) restricted to agents a and b, where each task in @S must be performed by one of
    % them and no other task is performed. Returns which tasks of @S go to agent a (@to_a) and the
    % solutions of both agents.
    ag_a = ev.mrta.agent{a};
    ag_b = ev.mrta.agent{b};
    P_a = ag_a.miqp.primal;
    P_b = ag_b.miqp.primal;
    n_a = ag_a.dv.nz;
    n_b = ag_b.dv.nz;
    th_a = ag_a.dv.n_e + ag_a.dv.n; % Offset of the task variables theta of agent a
    th_b = n_a + ag_b.dv.n_e + ag_b.dv.n; % Offset of the task variables theta of agent b (in the pair vector)

    % Bounds: theta = 0 for the tasks not in S
    lbz = [P_a.lbz; P_b.lbz];
    ubz = [P_a.ubz; P_b.ubz];
    ubz(th_a + (1:ag_a.dv.n_t)) = 0;
    ubz(th_b + (1:ag_b.dv.n_t)) = 0;
    loc_a = ev.idx_loc(S, a);
    loc_b = ev.idx_loc(S, b);
    ubz(th_a + loc_a(loc_a > 0)) = 1;
    ubz(th_b + loc_b(loc_b > 0)) = 1;

    % Coupling constraints (3) for the tasks in S
    n_S = length(S);
    rows = [find(loc_a > 0); find(loc_b > 0)];
    cols = [th_a + loc_a(loc_a > 0); th_b + loc_b(loc_b > 0)];
    C_S = sparse(rows, cols, 1, n_S, n_a + n_b);

    c = [ev.c_a{a}; ev.c_a{b}];
    A = blkdiag(sparse(P_a.A), sparse(P_b.A));
    C = [blkdiag(sparse(P_a.C), sparse(P_b.C)); C_S];
    miqp = MIQPholder([], c, A, [P_a.lb; P_b.lb], [P_a.ub; P_b.ub], C, [P_a.b; P_b.b; ones(n_S, 1)], ...
        [], [], [], lbz, ubz, [P_a.idx_int(:); n_a + P_b.idx_int(:)]);
    [x, V, info] = miqp.solve(ev.solver_args{:});
    ev.record_time(a, info); % Solved by agent a

    ok = info.feasible;
    to_a = false(n_S, 1);
    z_a = [];
    z_b = [];
    if ok
        z_a = ag_a.round_z(x(1:n_a));
        z_b = ag_b.round_z(x(n_a + (1:n_b)));
        to_a(loc_a > 0) = z_a(th_a + loc_a(loc_a > 0)) > 0.5;
        V = ev.c_a{a}'*z_a + ev.c_a{b}'*z_b;
    end
end

function z = selected_route_z(ag, route, S, ev, i)
    % Solution of agent @ag (index @i) that performs the tasks @S along the selected route @route, up to
    % the location of its last task, where it stops. Used to warm-start the solver.
    last = find(ismember(route, ev.task_v(S)), 1, 'last');
    if isempty(last); last = 1; end
    z = ag.route_to_z(route(1:last), ev.idx_loc(S, i));
end
