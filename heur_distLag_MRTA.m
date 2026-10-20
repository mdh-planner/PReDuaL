%% heur_distLag_MRTA - Heuristic dual Lagrangian method for MRTA problems (PReDuaL)
%
% Implements Algorithm 1 (heuristic subgradient method) and Algorithm 2 (heuristic update of the
% residual) of the article, followed by the primal recovery of Section VI (Algorithm 3):
%   B. Miloradovic, P. Krupa, M. Frasheri, C. Trubiani, "PReDuaL: A Primal-Recovering Residual-Driven
%   Dual Lagrangian Method for Multi-Robot Task Allocation".
%
% The coupling constraints (3) of the central MILP (9) are dualized, so that each agent solves its own
% MILP problem (13). Each agent s_i has its own dual variables lambda_i, which are updated with (16)
% using the residual returned by Algorithm 2. The preconditioning of Section VII is applied (cost
% normalization, sign inversion of the task costs and scaling of (10b) by 1/m).
% Set 'recovery' to false to run Algorithm 1 alone.
%
% INPUTS:
%   @mrta: Instance of MRTA. Its costs are normalized in place if 'normalize' is true
%   @a_loc: Cell array with the location of each agent
%   @t_loc: Cell array of cell arrays with the locations of the pending tasks of each task-type
%
% NAME/VALUE INPUTS (the defaults of the hyperparameters are the values of Section VII-A):
%   alpha: Step-size of Algorithm 1 (default 2*m*sum_k rho_{t_k})
%   alpha_rho: Multiplier (< 0) of the task costs in the inner MILPs (default -sum_{i,j} w_{i,j}/(2*sqrt(n)))
%   k_alpha: Number of iterations with constant step-size alpha (default 2)
%   k_h: Iteration from which Algorithm 2 is applied (default 5)
%   k_f: Over-feasible solutions are accepted for k > k_f (default 20)
%   k_max: Maximum number of iterations (default 50)
%   mu_up, mu_low: mu-bar and mu-underline of Algorithm 2 (default 1/m, 1)
%   gamma_up, gamma_low: gamma-bar and gamma-underline of Algorithm 2 (default 1, 1/m)
%   tol: Tolerance used to determine the sign of the elements of the residual (default 1e-6)
%   normalize: Apply the cost normalization of Section VII (default true)
%   recovery: Apply the primal recovery after Algorithm 1 (default true)
%   route_search: Apply the route search of the primal recovery (default true)
%   route_dp: In the primal recovery, try to compute the route of an agent for a set of tasks with
%             Held-Karp before solving its MILP problem (see RouteEvaluator.route_dp, default true)
%   recovery_rounds: Maximum number of rounds N_rec of the primal recovery (default 1). Without route
%                    search, a single round is used
%   column_generation: Start the primal recovery with column generation on the routes of Algorithm 1
%                      (default true), which also gives a lower bound of the optimal cost of (9) (info.LB).
%                      If the first route selection is within 1e-4 of the lower bound, the rest of the
%                      recovery is skipped
%   cg_stall: Stop the column generation if the value of the LP relaxation (19) did not decrease in the
%             last cg_stall rounds (k_cg). 0 to disable (default 3)
%   cg_max_time: Time limit [s] of the pricing MILP problems (20) of the column generation (default 2; NaN:
%                the time limit 'inner_max_time' of Algorithm 1). If it is reached, the best route found is
%                used and the lower bound uses the bound of the solver. When the column generation would
%                stop, it continues with the time limit 'inner_max_time' until it stops again
%   cg_max_iter: Maximum number of column generation rounds N_cg (default 20)
%   cg_smoothing: Dual smoothing parameter zeta of the column generation, in [0, 1) (default 0.5)
%   cg_tol: Relative tolerance epsilon_cg used to stop the column generation (default 1e-4)
%   tight: Add the constraints (7) to the MILP problems of the agents (default true). They do not change
%          the optimal routes, but make the LP relaxation much tighter (see Agent.gen_valid_ineq_constraints)
%   prize_bounds: At each iteration, fix to 0 the task variables with nonnegative cost and the stops at
%                 locations without tasks of negative cost (Section VII, default true)
%   warm_start: Warm-start the solver of each agent with its solution of the previous iteration (default true)
%   pool_stall: With the primal recovery, Algorithm 1 stops if no agent computed a new route in the last
%               pool_stall iterations (k_s in Section V). 0 to disable (default 3)
%   solver: Structure with name/value options for MIQPholder.solve (e.g., struct('solver', "intlinprog"))
%   inner_max_time: Time limit [s] of the MILP problems (13) of the agents in Algorithm 1 (default 10; inf:
%                   none). If it is reached, the best solution found by the solver is used (or the warm
%                   start, or staying at the initial location, if the solver found none). The re-routing
%                   problems of the primal recovery do not use this limit
%   parallel: Solve the MILP problems of the agents in parallel with parfor (default false). Requires the
%             Parallel Computing Toolbox
%   verbose: Print information at each iteration (default false)
%
% OUTPUTS:
%   @z_out: Solution z = [z_1; ...; z_m] (Algorithm 1 followed by the primal recovery)
%   @V_out: Cost (8) of @z_out (with the normalized costs)
%   @info: Structure with additional information:
%       .e_flag:  1 feasible solution (always the case if the primal recovery finds a solution)
%                 2 over-feasible solution found for k > k_f
%                 0 k_max (or pool_stall) reached: the best over-feasible solution found is returned
%                -1 k_max (or pool_stall) reached without an (over-)feasible solution: last iterate
%       .feasible: true if e_flag == 1
%       .e_flag_alg1, .V_alg1: Exit flag and cost of the solution of Algorithm 1 (before the recovery)
%       .recovery: Information of the primal recovery (see primal_recovery_MRTA)
%       .LB: Lower bound of the optimal cost of (9) from the column generation (-Inf if not computed)
%       .gap: Relative optimality gap (V_out - LB)/LB of @z_out (NaN if not computed)
%       .k: Number of iterations of Algorithm 1
%       .num_reused: Number of MILP problems that were not solved again because the costs and bounds
%                    of the agent had not changed
%       .num_inner_tl, .num_inner_nosol: Number of MILP problems stopped by 'inner_max_time' with and
%                                        without a solution
%       .stop_pool_stall: true if Algorithm 1 was stopped by 'pool_stall'
%       .run_time_par_total: Computation time of Section VII-A, assuming that the MILP problems of the
%                            agents are solved in parallel: overhead + sum_k max_i (time inner MILP i).
%                            Includes the primal recovery (computed in the same way)
%       .run_time_par_total_alg1: Computation time of Algorithm 1
%       .run_time_par: sum_k max_i (time inner MILP i)
%       .run_time_overhead: Time of Algorithm 1 except Steps 1 and 5
%       .run_time_build: Time of Step 1 (and of starting the parallel pool)
%       .run_time_all: Measured wall-clock time of Algorithm 1
%       .run_time_wall_total: Measured wall-clock time of Algorithm 1 plus the primal recovery
%       .parallel: true if the MILP problems were solved in parallel
%       .run_time_inner: m x k matrix with the solver time of each inner MILP
%       .time_best, .time_best_wall: Computation time (as .run_time_par_total) and wall-clock time at which
%                                    @z_out was found (NaN if no solution)
%       .time_best_alg1, .time_best_alg1_wall: The same for the solution of Algorithm 1
%       .best_hist: [time, wall-clock time, cost] of each new best feasible solution of (9) (Algorithm 1, then
%                   primal recovery)
%       .lambda: Dual variables (column i is lambda_i) of @z_out
%       .r: Residual (15) of @z_out
%       .hist: Residual (.r) and cost (.V) at each iteration
%       .hyper: Hyperparameters used
%       .c_max: Cost normalization factor (1 if 'normalize' is false)
%

function [z_out, V_out, info] = heur_distLag_MRTA(mrta, a_loc, t_loc, opt)
    arguments
        mrta MRTA
        a_loc
        t_loc
        opt.alpha (1, 1) double = NaN
        opt.alpha_rho (1, 1) double = NaN
        opt.k_alpha (1, 1) double {mustBeInteger, mustBeNonnegative} = 2
        opt.k_h (1, 1) double {mustBeInteger, mustBeNonnegative} = 5
        opt.k_f (1, 1) double {mustBeInteger, mustBeNonnegative} = 20
        opt.k_max (1, 1) double {mustBeInteger, mustBePositive} = 50
        opt.mu_up (1, 1) double = NaN
        opt.mu_low (1, 1) double = 1
        opt.gamma_up (1, 1) double = 1
        opt.gamma_low (1, 1) double = NaN
        opt.tol (1, 1) double {mustBePositive} = 1e-6
        opt.normalize (1, 1) {mustBeNumericOrLogical} = true
        opt.recovery (1, 1) {mustBeNumericOrLogical} = true
        opt.route_search (1, 1) {mustBeNumericOrLogical} = true
        opt.route_dp (1, 1) {mustBeNumericOrLogical} = true
        opt.recovery_rounds (1, 1) double {mustBeInteger, mustBePositive} = 1
        opt.column_generation (1, 1) {mustBeNumericOrLogical} = true
        opt.cg_max_iter (1, 1) double {mustBeInteger, mustBeNonnegative} = 20
        opt.cg_smoothing (1, 1) double {mustBeNonnegative} = 0.5
        opt.cg_tol (1, 1) double {mustBePositive} = 1e-4
        opt.cg_stall (1, 1) double {mustBeInteger, mustBeNonnegative} = 3
        opt.cg_max_time (1, 1) double = 2
        opt.tight (1, 1) {mustBeNumericOrLogical} = true
        opt.prize_bounds (1, 1) {mustBeNumericOrLogical} = true
        opt.warm_start (1, 1) {mustBeNumericOrLogical} = true
        opt.pool_stall (1, 1) double {mustBeInteger, mustBeNonnegative} = 3
        opt.solver struct = struct()
        opt.inner_max_time (1, 1) double {mustBePositive} = 10
        opt.parallel (1, 1) {mustBeNumericOrLogical} = false
        opt.verbose (1, 1) {mustBeNumericOrLogical} = false
    end

    m = mrta.m;
    start_build = tic;

    %% Step 1: Construct the ingredients of problem (10)

    % Cost normalization (Section VII)
    info.c_max = 1;
    if opt.normalize
        info.c_max = mrta.normalize_costs();
    end

    % MILP problem (13) of each agent: sets Z_i and costs q_i
    mrta.build_distMILP(a_loc, t_loc, 'condense', true, 'tight', opt.tight);

    % Coupling constraints (3), written as sum_i A_i*z_i = b in (10b)
    [A, b, idx_A] = mrta.gen_single_task_eq_constraints(t_loc);
    A = A/m; % Preconditioning of (10b): all non-zero elements are equal to 1/m
    b = b/m;
    m_z = size(A, 1); % Number of coupling constraints (pending tasks)

    % Default values of the hyperparameters (Section VII-A)
    if isnan(opt.alpha); opt.alpha = 2*m*sum(mrta.param.cost_task); end
    if isnan(opt.alpha_rho); opt.alpha_rho = -sum(mrta.Gv.Edges.Weight)/(2*sqrt(mrta.n)); end
    if isnan(opt.mu_up); opt.mu_up = 1/m; end
    if isnan(opt.gamma_low); opt.gamma_low = 1/m; end
    if opt.alpha_rho >= 0
        warning("heur_distLag_MRTA:alpha_rho", "alpha_rho should be negative (see Section V-A)");
    end
    solver_args = namedargs2cell(opt.solver);
    % If intlinprog is used and the MOSEK toolbox shadows it, remove MOSEK from the path until we return
    optim_guard = MIQPholder.optim_path_guard(solver_args); %#ok<NASGU>
    % Solver options of the MILP problems of the agents in Algorithm 1 (time limit 'inner_max_time')
    alg1_solver = opt.solver;
    if isfinite(opt.inner_max_time)
        if isfield(alg1_solver, 'max_time')
            alg1_solver.max_time = min(alg1_solver.max_time, opt.inner_max_time);
        else
            alg1_solver.max_time = opt.inner_max_time;
        end
    end
    alg1_args = namedargs2cell(alg1_solver);

    % Ingredients of each agent
    A_a = cell(1, m); % A_i
    c_a = cell(1, m); % Costs of (8)
    q_a = cell(1, m); % Costs of (8) with the task costs multiplied by alpha_rho (Section V-A)
    for i = 1:m
        ag = mrta.agent{i};
        A_a{i} = A(:, ag.dv.idx + (1:ag.dv.nz));
        c_a{i} = ag.miqp.primal.c;
        q_a{i} = c_a{i};
        idx_theta = ag.dv.n_e + ag.dv.n + (1:ag.dv.n_t);
        q_a{i}(idx_theta) = opt.alpha_rho*q_a{i}(idx_theta);
    end
    c = cat(1, c_a{:});

    % Solver of the MILP problems of the agents (sequential or parallel)
    agent_solver = AgentSolver(mrta, alg1_args, opt.parallel);
    info.parallel = agent_solver.parallel;
    info.run_time_build = toc(start_build);

    %% Steps 2-18: heuristic subgradient method
    start_loop = tic;
    lambda = zeros(m_z, m); % Step 2: lambda(:, i) is lambda_i
    Gam = inf(m_z, m); % Gam(j, i) = Gamma_i(t_j) at the current iteration (Inf if s_i does not do t_j)
    Gam_best = inf(m_z, m); % Smallest Gamma_i(t_j) obtained in past iterations
    z_a = cell(1, m);
    t_inner = zeros(m, opt.k_max); % Solver time of each inner MILP
    t_batch = zeros(1, opt.k_max); % Wall time of Steps 4-6 at each iteration
    t_iter = zeros(1, opt.k_max); % Wall time from the start of the loop to the solution of each iteration
    info.hist.r = zeros(m_z, opt.k_max);
    info.hist.V = zeros(1, opt.k_max);
    z_best = [];
    k_best = -1; % Iteration of z_best
    V_best = inf;
    done = false;
    pool = cell(1, m); % Routes computed by each agent (route pool of the primal recovery)
    for i = 1:m
        pool{i} = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    stall = 0; % Consecutive iterations without new routes (see 'pool_stall')
    info.stop_pool_stall = false;
    c_prev = cell(1, m); % Costs and bounds of the problem of each agent at the previous iteration
    lbz_prev = cell(1, m);
    ubz_prev = cell(1, m);
    z_prev = cell(1, m); % Solution and information of the problem of each agent at the previous iteration
    info_prev = cell(1, m);
    info.num_reused = 0;
    info.num_inner_tl = 0;
    info.num_inner_nosol = 0;

    for k = 0:(opt.k_max - 1)

        % Steps 4-6: each agent solves its MILP problem (13)
        c_k = cell(1, m);
        lbz_k = cell(1, m); % Bounds of the variables (default bounds of the agent if empty)
        ubz_k = cell(1, m);
        ws_k = cell(1, m); % Initial solutions for the solver (none if empty)
        for i = 1:m
            c_k{i} = q_a{i} + A_a{i}'*lambda(:, i);
            ag = mrta.agent{i};
            if opt.prize_bounds
                [lbz_k{i}, ubz_k{i}] = ag.prize_bounds(c_k{i}, ag.miqp.primal.lbz, ag.miqp.primal.ubz);
            end
            if opt.warm_start && k > 0
                ws_k{i} = ag.warm_start_z(z_a{i}, c_k{i});
            end
        end
        % An agent whose problem did not change since the previous iteration keeps its previous solution
        solve_k = true(1, m);
        if k > 0
            for i = 1:m
                solve_k(i) = ~(isequal(c_k{i}, c_prev{i}) && isequal(lbz_k{i}, lbz_prev{i}) && isequal(ubz_k{i}, ubz_prev{i}));
            end
        end
        z_k_list = z_prev;
        info_k_list = info_prev;
        for i = find(~solve_k)
            info_k_list{i}.run_time = 0;
        end
        [z_k_list(solve_k), ~, info_k_list(solve_k), t_batch(k+1)] = ...
            agent_solver.solve(find(solve_k), c_k(solve_k), lbz_k(solve_k), ubz_k(solve_k), ws_k(solve_k));
        info.num_reused = info.num_reused + nnz(~solve_k);
        z_prev(solve_k) = z_k_list(solve_k);
        info_prev(solve_k) = info_k_list(solve_k);
        c_prev = c_k;
        lbz_prev = lbz_k;
        ubz_prev = ubz_k;

        num_new = 0; % Number of new routes in the pool
        for i = 1:m
            ag = mrta.agent{i};
            info_i = info_k_list{i};
            if solve_k(i) && ~info_i.feasible && isfinite(opt.inner_max_time)
                % Time limit reached without a solution: use the warm start, or stay at the initial location
                if ~isempty(ws_k{i})
                    z_k_list{i} = ws_k{i};
                else
                    z_k_list{i} = ag.route_to_z(ag.loc, []);
                end
                info_i.feasible = true;
                info_k_list{i} = info_i;
                z_prev{i} = z_k_list{i};
                info_prev{i} = info_i;
                info.num_inner_nosol = info.num_inner_nosol + 1;
            elseif solve_k(i) && info_i.e_flag == 2
                info.num_inner_tl = info.num_inner_tl + 1;
            end
            if ~info_i.feasible
                details = "";
                if isfield(info_i, 'mosek')
                    details = sprintf(" (MOSEK response code %d: %s, problem status: %s, solution status: %s, max. violation: %g)", ...
                        info_i.mosek_e_flag, info_i.mosek.rcodestr, info_i.mosek.prosta, info_i.mosek.solsta, info_i.mosek.violation);
                elseif isfield(info_i, 'intlinprog')
                    details = sprintf(" (intlinprog exit flag %d)", info_i.intlinprog.exitflag);
                end
                error("heur_distLag_MRTA:inner_MILP", "Inner MILP of agent %s could not be solved at iteration %d%s", ag.name, k, details);
            end
            t_inner(i, k+1) = info_i.run_time;
            z_a{i} = ag.round_z(z_k_list{i});

            % Time Gamma_i(t_j) at which the agent completes each task t_j (Section V-B)
            Gam_i = ag.task_times(z_a{i});
            Gam(:, i) = Gam_i(sub2ind(size(Gam_i), idx_A.v(:), idx_A.t(:)));

            % Add the route of the agent to the pool
            route = ag.get_path(z_a{i});
            key = sprintf('%d,', route);
            if ~isKey(pool{i}, key)
                pool{i}(key) = route;
                num_new = num_new + 1;
            end
        end
        Gam(Gam == 0) = inf; % Tasks that the agent cannot do
        Gam_best = min(Gam_best, Gam);

        % Step 7: residual (15)
        z_k = cat(1, z_a{:});
        r = A*z_k - b;
        V_k = c'*z_k;
        info.hist.r(:, k+1) = r;
        info.hist.V(k+1) = V_k;
        lambda_k = lambda;

        % Feasible, over-feasible or infeasible
        is_feas = all(abs(r) <= opt.tol);
        is_over = ~is_feas && all(r >= -opt.tol);
        if is_over && V_k < V_best
            z_best = z_k;
            V_best = V_k;
            lambda_best = lambda_k;
            r_best = r;
            k_best = k;
        end
        t_iter(k+1) = toc(start_loop);

        if opt.verbose
            fprintf("k = %3d | J = %10.4f | tasks done once: %d, more than once: %d, not done: %d\n", ...
                k, V_k, nnz(abs(r) <= opt.tol), nnz(r > opt.tol), nnz(r < -opt.tol));
        end

        % Termination conditions
        if is_feas
            info.e_flag = 1;
            done = true;
            break;
        elseif is_over && k > opt.k_f
            info.e_flag = 2;
            done = true;
            break;
        end
        if num_new == 0
            stall = stall + 1;
        else
            stall = 0;
        end
        if opt.recovery && opt.pool_stall > 0 && stall >= opt.pool_stall
            info.stop_pool_stall = true;
            break; % The best over-feasible solution (or the last iterate) is returned
        end

        % Steps 8-14: heuristic update of the residual (Algorithm 2)
        if k >= opt.k_h
            R = heuristic_residual(r, Gam, Gam_best, opt);
        else
            R = repmat(r, 1, m);
        end

        % Steps 15-16: step-size
        if k > opt.k_alpha
            alpha_k = opt.alpha/(k - opt.k_alpha);
        else
            alpha_k = opt.alpha;
        end

        % Step 17: update of the dual variables (16)
        lambda = lambda + alpha_k*R;

    end

    % Stopped without a feasible solution: return the best over-feasible one (if any)
    k_sol = k; % Iteration of the returned solution
    if ~done
        if ~isempty(z_best)
            info.e_flag = 0;
            z_k = z_best;
            lambda_k = lambda_best;
            r = r_best;
            k_sol = k_best;
        else
            info.e_flag = -1;
        end
    end
    info.run_time_all = toc(start_loop);

    %% Solution of Algorithm 1
    num_iters = k + 1;
    z_out = mrta.round_res(z_k);
    V_out = c'*z_out;

    info.feasible = info.e_flag == 1;
    info.k = num_iters;
    info.lambda = lambda_k;
    info.r = r;
    info.idx_A = idx_A;
    info.hist.r = info.hist.r(:, 1:num_iters);
    info.hist.V = info.hist.V(1:num_iters);
    info.hyper = rmfield(opt, {'solver', 'verbose'});
    info.solver = opt.solver;

    % Computation times (Section VII-A)
    info.run_time_inner = t_inner(:, 1:num_iters);
    info.run_time_par = sum(max(info.run_time_inner, [], 1));
    info.run_time_overhead = max(info.run_time_all - sum(t_batch(1:num_iters)), 0);
    info.run_time_par_total = info.run_time_overhead + info.run_time_par;
    info.run_time_seq = sum(info.run_time_inner(:));
    info.run_time_par_total_alg1 = info.run_time_par_total;
    info.e_flag_alg1 = info.e_flag;
    info.V_alg1 = V_out;

    % Time at which the solution of Algorithm 1 was found (computed as info.run_time_par_total)
    info.time_best_alg1 = NaN;
    info.time_best_alg1_wall = NaN;
    if info.e_flag >= 0
        info.time_best_alg1 = sum(max(t_inner(:, 1:k_sol+1), [], 1)) + max(t_iter(k_sol+1) - sum(t_batch(1:k_sol+1)), 0);
        info.time_best_alg1_wall = t_iter(k_sol+1);
    end
    info.best_hist = zeros(0, 3);
    if info.e_flag == 1
        info.best_hist(1, :) = [info.time_best_alg1, info.time_best_alg1_wall, V_out];
    end
    info.LB = -inf;
    info.gap = NaN;

    %% Primal recovery (Section VI)
    info.run_time_wall_total = info.run_time_all;
    if opt.recovery && m_z > 0
        rec_opt.rounds = 1;
        if opt.route_search
            rec_opt.rounds = opt.recovery_rounds;
        end
        rec_opt.route_search = opt.route_search;
        rec_opt.route_dp = opt.route_dp;
        rec_opt.column_generation = opt.column_generation;
        rec_opt.cg_max_iter = opt.cg_max_iter;
        rec_opt.cg_smoothing = opt.cg_smoothing;
        rec_opt.cg_tol = opt.cg_tol;
        rec_opt.cg_stall = opt.cg_stall;
        % Solver options of the pricing problems: time limit of Algorithm 1 unless cg_max_time is given
        cg_solver = alg1_solver;
        if ~isnan(opt.cg_max_time)
            cg_solver.max_time = opt.cg_max_time;
        end
        rec_opt.cg_solver_args = namedargs2cell(cg_solver);
        if ~isnan(opt.cg_max_time)
            rec_opt.cg_solver_args_final = alg1_args; % Last rounds with the time limit of Algorithm 1
        end
        rec_opt.solver_args = solver_args;
        agent_solver.solver_args = solver_args; % No 'inner_max_time' in the primal recovery
        rec_opt.agent_solver = agent_solver;
        [z_rec, V_rec, info.recovery] = primal_recovery_MRTA(mrta, idx_A, pool, Gam, rec_opt);
        if ~isempty(z_rec) && (info.e_flag ~= 1 || V_rec <= V_out)
            z_out = z_rec;
            V_out = V_rec;
            info.e_flag = 1;
            info.feasible = true;
            info.r = A*z_out - b;
        end
        info.LB = info.recovery.LB;
        if isfinite(info.LB)
            info.gap = (V_out - info.LB)/max(abs(info.LB), 1e-12);
        end
        % New best solutions of the recovery (times from the start of Algorithm 1)
        for q = 1:size(info.recovery.best_hist, 1)
            row = info.recovery.best_hist(q, :) + [info.run_time_par_total_alg1, info.run_time_all, 0];
            if isempty(info.best_hist) || row(3) < info.best_hist(end, 3) - 1e-9
                info.best_hist(end+1, :) = row;
            end
        end
        info.run_time_par_total = info.run_time_par_total + info.recovery.run_time_par_total;
        info.run_time_wall_total = info.run_time_wall_total + info.recovery.run_time_all;
        info.run_time_seq = info.run_time_seq + info.recovery.run_time_seq;
        if opt.verbose
            fprintf("Primal recovery: cost %.4f -> %.4f (%d MILPs solved by the agents, %d route selection problems, %d routes)\n", ...
                info.V_alg1, V_out, info.recovery.num_solves, info.recovery.num_selections, info.recovery.pool_size);
            if isfinite(info.LB)
                fprintf("Lower bound: %.4f, certified optimality gap: %.2f%% (%d column generation rounds)\n", ...
                    info.LB, 100*info.gap, info.recovery.cg.iters);
            end
        end
    end

    % Time at which @z_out was found: first new best feasible solution with its cost (or the solution of
    % Algorithm 1, if it is returned without being feasible)
    info.time_best = info.time_best_alg1;
    info.time_best_wall = info.time_best_alg1_wall;
    if info.feasible
        q = find(info.best_hist(:, 3) <= V_out + 1e-9*max(1, abs(V_out)), 1);
        if ~isempty(q)
            info.time_best = info.best_hist(q, 1);
            info.time_best_wall = info.best_hist(q, 2);
        end
    end

end

function R = heuristic_residual(r, Gam, Gam_best, opt)
    % Algorithm 2: heuristic update of the residual.
    % Returns the matrix R whose column i is hat{r}_i = r + Delta r_i.
    %   @r: Residual (15)
    %   @Gam: Gam(j, i) = Gamma_i(t_j) for the current solutions (Inf if s_i does not do t_j)
    %   @Gam_best: Smallest Gamma_i(t_j) obtained in past iterations (Inf if never done)
    %   @opt: Structure with the hyperparameters
    [m_z, m] = size(Gam);
    dR = zeros(m_z, m);

    for j = 1:m_z
        does = isfinite(Gam(j, :)); % Agents whose current solution performs task t_j

        if abs(r(j)) <= opt.tol
            % Steps 4-6: task performed by exactly one agent, the residual is not modified
            continue

        elseif r(j) > 0
            % Steps 7-16: task performed by several agents.
            % s_l is the agent that completes t_j in the least time (lowest index in case of ties)
            [G_min, l] = min(Gam(j, :));
            G_max = max(Gam(j, does));
            dT = zeros(1, m);
            dT(does) = (Gam(j, does) - G_min)/G_max; % Delta T_{i,j}^[k], see (17)
            dR(j, :) = opt.mu_up*dT; % Disincentivize the other agents according to (17)
            dR(j, l) = -opt.mu_low*r(j); % Incentivize s_l

        else
            % Steps 17-23: task not performed by any agent.
            % s_l is the agent that has completed t_j in the least time in past iterations.
            % If no agent has ever performed t_j, then Delta r_i(j) = 0 for all agents.
            [G_min, l] = min(Gam_best(j, :));
            if isfinite(G_min)
                dR(j, :) = -opt.gamma_low*r(j); % Disincentivize the other agents
                dR(j, l) = opt.gamma_up*r(j); % Incentivize s_l
            end
        end
    end

    R = r + dR;
end
