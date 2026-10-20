%% primal_recovery_MRTA - Primal recovery of PReDuaL (Section VI and Algorithm 3 of the article)
%
% Builds a feasible solution of the central MILP (9) from the routes computed by the agents during
% Algorithm 1 (route pool).
%
% Column generation (opt.column_generation): the LP relaxation (19) of the route selection problem (18)
% is solved over the pool, and the agents price new routes with their MILP problems (20) using its dual
% variables, which also gives the Lagrangian lower bound (21) of the optimal cost of (9). If the route
% selection is then within opt.gap_tol of the lower bound, its routes are returned directly.
%
% Then, each round of the recovery does the following:
%   1. Route selection: solve problem (18), which picks one route of the pool per agent and assigns
%      each pending task to an agent whose route visits its location (see RoutePool.selection_problem).
%   2. Route search (opt.route_search): generate new routes from the selected ones without solving any
%      MILP problem (routes (a)-(d) of Section VI-C, see RouteGraph), add them to the pool and solve (18)
%      again, while its cost decreases.
%   3. Re-routing: each agent computes its optimal route for the tasks assigned to it.
% The routes computed in step 3 are added to the pool, and the rounds are repeated while the cost
% decreases (at most opt.rounds times). If (18) is infeasible, the last iterate of Algorithm 1 is
% repaired instead: a task done by several agents is kept by the one that completes it first, and a
% task done by none is inserted where it is cheapest.
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
%       .column_generation: (optional) Apply the column generation before the rounds (default false)
%       .cg_max_iter: Maximum number of column generation rounds
%       .cg_smoothing: Dual smoothing parameter in [0, 1) (0 for no smoothing)
%       .cg_tol: Relative tolerance used to stop the column generation
%       .cg_stall: (optional) Stop the column generation if the value of the LP relaxation (19) did not
%                  decrease (relative tolerance cg_tol) in the last cg_stall rounds. 0 to disable (default 0)
%       .cg_solver_args: (optional) Options for MIQPholder.solve of the pricing problems of the agents, e.g.,
%                        with a time limit (default: solver_args)
%       .cg_solver_args_final: (optional) Options of the pricing problems used once the column generation
%                              would stop with cg_solver_args (e.g., with a longer time limit): the column
%                              generation then continues with them until it stops again (default: none)
%       .gap_tol: (optional) Relative tolerance of the test after the column generation (default 1e-4)
%       .solver_args: Cell array with name/value options for MIQPholder.solve
%       .agent_solver: (optional) AgentSolver for the MILP problems of the agents (sequential if not given)
%
% OUTPUTS:
%   @z: Feasible solution of (9) (empty if none was found)
%   @V: Cost (8) of @z (Inf if no feasible solution was found)
%   @rec: Structure with information:
%       .LB: Lower bound of the optimal cost of (9) (-Inf if not computed)
%       .gap: Relative optimality gap (V - LB)/LB
%       .V_select, .V_search, .V_reroute: Cost after each step of each round
%       .used_repair: true if the last iterate of Algorithm 1 was repaired
%       .certified: true if, after the column generation, the first route selection was within opt.gap_tol of
%                   the lower bound (and was returned)
%       .cg: Information of the column generation
%       .num_solves: Number of MILP problems solved by the agents
%       .num_dp: Number of routes computed by RouteEvaluator.route_dp
%       .num_selections: Number of route selection problems solved
%       .pool_size: Number of routes in the pool at the end
%       .run_time_search: Time spent generating the routes of the route search
%       .run_time_par_total: Computation time, as in Section VI-D
%       .best_hist: [time, wall-clock time, cost] of each new best solution, from the start of the recovery
%                   (time computed as .run_time_par_total up to that point)
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
    rec.used_repair = false;
    rec.run_time_coord = 0; % Solver time of the route selection problems (done by the coordinator)
    rec.wall_time_coord = 0;
    rec.run_time_search = 0;
    rec.num_selections = 0;
    rec.best_hist = zeros(0, 3);
    rec.tic0 = start_time;

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

    if ~isfield(opt, 'column_generation'); opt.column_generation = false; end
    if ~isfield(opt, 'gap_tol'); opt.gap_tol = 1e-4; end
    rec.LB = -inf;
    rec.cg = struct('iters', 0, 'num_new', 0, 'V_lp', NaN, 'converged', false, 'stalled', false, ...
                    'hist_LB', zeros(1, 0), 'hist_V_lp', zeros(1, 0), 'final_rounds', 0, 'run_time', 0, ...
                    'time_lp', 0, 'time_pricing', 0);
    rec.certified = false;

    %% Column generation on the pool of Algorithm 1
    if opt.column_generation
        [LB_cg, cg, rec] = column_generation(rp, ev, opt, rec);
        rec.LB = LB_cg;
        rec.cg = cg;
        % The best combination of routes may already be optimal up to gap_tol
        [sel0, ok0, rec] = select_routes(rp, ev, rec);
        if ok0
            rec.V_select(end+1) = sel0.V;
            rec.certified = sel0.V - rec.LB <= opt.gap_tol*max(1, abs(sel0.V));
        end
    end

    %% Rounds of route selection, route search and re-routing
    if rec.certified
        assign = sel0.assign;
        V = sel0.V;
        rec.best_hist(end+1, :) = [paper_time(rec, ev), toc(rec.tic0), V];
    else
        [assign, V, rec] = improve(mrta, rp, graph, Gam, ev, opt, rec);
    end

    %% Construct the solution
    if isempty(assign)
        z = [];
        V = inf;
    else
        z = zeros(mrta.nz, 1);
        V = 0;
        for i = 1:m
            ag = mrta.agent{i};
            S_i = find(assign == i);
            z_i = [];
            if rec.certified
                % Selected route (no MILP problem to solve)
                z_i = ag.route_to_z(rp.paths{i}{sel0.route(i)}, ev.idx_loc(S_i, i));
                if ~ag.miqp.is_feasible(z_i); z_i = []; end
            end
            if isempty(z_i)
                [~, z_i] = ev.cost(i, S_i);
            end
            z(ag.dv.idx + (1:ag.dv.nz)) = z_i;
            V = V + ev.c_a{i}'*z_i;
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
    rec = rmfield(rec, 'tic0');

end

function t = paper_time(rec, ev)
    % Computation time of the recovery so far, as rec.run_time_par_total (Section VI-D)
    wall = toc(rec.tic0);
    t = rec.run_time_coord + ev.run_time_par + max([ev.round_times, 0]) + max(wall - rec.wall_time_coord - ev.run_time_wall, 0);
end

function [assign_best, V_best, rec] = improve(mrta, rp, graph, Gam, ev, opt, rec)
    % Rounds of the primal recovery, repeated while the cost decreases. Returns the best assignment
    % @assign_best and its cost @V_best (empty and Inf if none was found).
    m = mrta.m;
    assign_best = [];
    V_best = inf;
    searched = containers.Map('KeyType', 'char', 'ValueType', 'logical'); % See route_search
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
        elseif r_idx == 1
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

        % Add the new routes to the pool
        num_new = add_routes_to_pool(rp, ev);

        if V_round < V_best - 1e-9
            V_best = V_round;
            assign_best = assign;
            rec.best_hist(end+1, :) = [paper_time(rec, ev), toc(rec.tic0), V_best];
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

function [LB, cg, rec] = column_generation(rp, ev, opt, rec)
    % Column generation on the LP relaxation of the route selection problem. In each round, the LP
    % relaxation is solved over the pool, and the agents price new routes with its dual variables (each
    % agent is warm-started with its route of the LP with the smallest reduced cost).
    LB = -inf;
    duals_best = [];
    cg = struct('iters', 0, 'num_new', 0, 'V_lp', NaN, 'converged', false, 'stalled', false, ...
                'hist_LB', zeros(1, 0), 'hist_V_lp', zeros(1, 0), 'final_rounds', 0, 'run_time', 0, ...
                'time_lp', 0, 'time_pricing', 0);
    if ~isfield(opt, 'cg_stall'); opt.cg_stall = 0; end
    tol_lp = NaN; % Set with the value of the first LP
    time_start = rec.run_time_coord + ev.run_time_par;
    % Solver options of the pricing problems (e.g., with a time limit), restored at the end
    solver_args_prev = ev.solver.solver_args;
    if isfield(opt, 'cg_solver_args') && ~isempty(opt.cg_solver_args)
        ev.solver.solver_args = opt.cg_solver_args;
    end
    args_final = {};
    if isfield(opt, 'cg_solver_args_final'); args_final = opt.cg_solver_args_final; end
    stall_from = 1; % First LP of the stall test

    for it = 1:opt.cg_max_iter
        cg.iters = it;

        % LP relaxation of the route selection problem and its dual variables
        [duals, V_lp, ok, info_lp, lp] = selection_lp(rp, ev);
        rec.run_time_coord = rec.run_time_coord + info_lp.run_time;
        rec.wall_time_coord = rec.wall_time_coord + info_lp.wall_time;
        cg.time_lp = cg.time_lp + info_lp.run_time;
        if ~ok; break; end
        cg.V_lp = V_lp;
        cg.hist_V_lp(end+1) = V_lp;
        if isnan(tol_lp); tol_lp = opt.cg_tol*max(1, abs(V_lp)); end

        % Pricing problems of the agents: first with the smoothed duals and, if they give no new route,
        % with the duals of the LP
        candidates = {duals};
        if ~isempty(duals_best) && opt.cg_smoothing > 0
            candidates = {opt.cg_smoothing*duals_best + (1 - opt.cg_smoothing)*duals, duals};
        end
        num_new = 0;
        t_par = ev.run_time_par;
        for q = 1:length(candidates)
            [L, valid, num_new] = price_all(candidates{q}, rp, ev, lp.best_route);
            if valid && L > LB
                LB = L;
                duals_best = candidates{q};
            end
            if num_new > 0; break; end
        end
        cg.time_pricing = cg.time_pricing + ev.run_time_par - t_par;
        cg.num_new = cg.num_new + num_new;
        cg.hist_LB(end+1) = LB;

        % Stopping criteria. Before stopping because no route improves the LP, the pricing problems are
        % solved with the final solver options (if given)
        if V_lp - LB <= tol_lp
            cg.converged = true; % LP relaxation solved: its value is the best Lagrangian bound
            break;
        end
        stalled = opt.cg_stall > 0 && length(cg.hist_V_lp) - stall_from >= opt.cg_stall && ...
                  cg.hist_V_lp(end - opt.cg_stall) - V_lp <= tol_lp;
        if num_new == 0 || stalled
            if ~isempty(args_final)
                ev.solver.solver_args = args_final;
                args_final = {};
                cg.final_rounds = it;
                stall_from = length(cg.hist_V_lp); % The stall test restarts with the final options
                continue
            end
            cg.converged = num_new == 0;
            cg.stalled = stalled && num_new > 0;
            break;
        end
    end
    ev.solver.solver_args = solver_args_prev;
    cg.run_time = rec.run_time_coord + ev.run_time_par - time_start;
end

function [L, valid, num_new] = price_all(duals, rp, ev, ws_routes)
    % Solves the pricing problem of each agent for the dual variables @duals, adds the new routes to the
    % pool and returns the Lagrangian lower bound @L (valid if the solver returned a finite lower bound
    % for all the problems). @ws_routes: (optional) route of each agent used to warm-start its solver
    if nargin < 4; ws_routes = {}; end
    [~, bound, z, info] = ev.price_agents(duals, ws_routes);
    L = sum(duals) + sum(bound);
    valid = all(isfinite(bound)); % The bounds of the solver are valid even if it stopped at a time limit
    num_new = 0;
    for i = 1:rp.m
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

function [duals, V, ok, info, lp] = selection_lp(rp, ev)
    % Solves the LP relaxation of the route selection problem and returns dual variables pi_j of the
    % constraints "task j is performed by exactly one agent" (one per pending task): pi_j = rho_{t_j} +
    % sigma_j, where sigma_j >= 0 is the dual variable of the covering constraint of task j (see
    % RoutePool.selection_problem). The pricing problems then use the task costs rho_{t_j} - pi_j = -sigma_j.
    % @lp: structure with the dual variables .sigma (covering) and .mu (one route per agent) and the route
    %      of each agent used by the LP with the smallest reduced cost (.best_route, to warm-start its
    %      pricing problem)
    [miqp, col_agent, col_route, offset] = rp.selection_problem(ev.task_rho);
    [u, V, info] = miqp.solve(ev.solver_args{:}, 'relax', true);
    V = V + offset;
    ok = info.feasible && isfield(info, 'y_in');
    duals = [];
    lp = struct('sigma', [], 'mu', [], 'best_route', {{}});
    if ok
        duals = ev.task_rho(:) + info.y_in(:);
        lp.sigma = max(info.y_in(:), 0);
        lp.mu = info.y_eq(:);
        lp.best_route = cell(1, rp.m);
        for i = 1:rp.m
            cols = find(col_agent == i & u(:) > 1e-6);
            best = inf;
            for c = reshape(cols, 1, [])
                rc_c = rp.E{i}(col_route(c)) - lp.mu(i) - sum(lp.sigma(rp.cov{i}{col_route(c)}));
                if rc_c < best
                    best = rc_c;
                    lp.best_route{i} = rp.paths{i}{col_route(c)};
                end
            end
        end
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

function z = selected_route_z(ag, route, S, ev, i)
    % Solution of agent @ag (index @i) that performs the tasks @S along the selected route @route, up to
    % the location of its last task, where it stops. Used to warm-start the solver.
    last = find(ismember(route, ev.task_v(S)), 1, 'last');
    if isempty(last); last = 1; end
    z = ag.route_to_z(route(1:last), ev.idx_loc(S, i));
end
