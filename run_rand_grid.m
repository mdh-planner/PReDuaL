%% run_rand_grid() - Generates and solves a random grid problem (Section VII-B of the article)
%
% The problem is solved with PReDuaL (heur_distLag_MRTA) and as the central MILP (9).
%
% INPUTS:
% @param: struct, containing the parameters for random generation (see gen_param_grid)
%   .row_lim: [int, int], min/max limits for number of rows of grid
%   .col_lim: [int, int], min/max limits for number of columns of grid
%   .m_lim: [int, int], min/max limits for number of agents
%   .k_lim: [int, int], min/max limits for number of task-types
%   .Evt_sparse_lim: [float, float] in [0, 1], sparsity of Evt
%   .Eat_sparse_lim: [float, float] in [0, 1], sparsity of Eat
%   .max_tv_lim: [int, int], min/max limits for max number of each task
%   .min_tv_lim: [int, int], min/max limits for min number of each task
% @seed: struct, used to set seeds for the random generation
%   .graph: int, seed for generation of random graph
%   .state: int, seed for assignment of tasks and locations
% @save_me: bool (optional), determines if the results are saved in a file
% @save_name: str (optional) no spaces, used for name generation of saved file
%
% NAME/VALUE INPUTS:
% 'save_to_results': bool, saves the mat files to results/
% 'save_add_time': bool, adds a timestamp to the name of the file
% 'run_heur': bool, determines if PReDuaL is run
% 'run_central': bool, determines if the central solver is run
% 'ignore_seed_central': array of state seeds for which the central solver is not run
% 'ignore_seed_heur': array of state seeds for which PReDuaL is not run
% 'heur': struct, name/value options for heur_distLag_MRTA
% 'solver': struct, name/value options for MIQPholder.solve (used for all MILPs), e.g. struct('solver', "intlinprog")
% 'central_max_time': maximum time [s] for the central MILP solver (default inf)
% 'central_tight': add the constraints (7) to the central MILP (9) (default false)
%
% OUTPUTS:
% @res: Structure containing results of the solver (timing results, num. iterations, exit flag, etc.)
% @test_info: Structure containing the information defining the test.
%             This information is saved to a .mat file in 'results/' if save_me = true
%

function [res, test_info] = run_rand_grid(param, seed, save_me, save_name, opt)
    arguments
        param struct
        seed struct
        save_me (1, 1) {mustBeNumericOrLogical} = true
        save_name {mustBeText} = ""
        opt.save_to_results = true
        opt.save_add_time = true
        opt.run_heur = true
        opt.run_central = true
        opt.ignore_seed_central = []
        opt.ignore_seed_heur = []
        opt.heur struct = struct()
        opt.solver struct = struct()
        opt.central_max_time (1, 1) double {mustBePositive} = inf
        opt.central_tight (1, 1) {mustBeNumericOrLogical} = false
    end
    res = struct();
    test_info = struct();

    % Check arguments
    if strcmp(save_name, ""); save_name = "test"; end
    run_heur = opt.run_heur && ~any(opt.ignore_seed_heur == seed.state);
    run_central = opt.run_central && ~any(opt.ignore_seed_central == seed.state);
    if ~run_heur && ~run_central
        fprintf(" Problem ignored! Graph seed: %d, State seed: %d\n\n", seed.graph, seed.state);
        return
    end

    %% Generate random grid (with the costs normalized as in Section VII-A)
    [prob, a_loc, t_loc, gen_info] = gen_rand_grid_instance(param, seed);
    n_rows = gen_info.n_rows; % Rows of grid
    n_cols = gen_info.n_cols; % Columns of grid
    num_tasks_tt = gen_info.num_tasks;
    prob.condense_Evt(t_loc); % Number of decision variables nz of the condensed problem (9)
    fprintf(" Problem with n = %d, m = %d, k = %d, num_tasks = %d, nz = %d.", prob.n, prob.m, prob.k, num_tasks_tt, prob.nz);
    fprintf(" Graph seed: %d, State seed: %d\n", seed.graph, seed.state);

    %% Save initial configuration
    save_name = sprintf("%s_gs_%s_ss_%s", save_name, num2str(seed.graph), num2str(seed.state));
    if opt.save_to_results
        if ~isfolder("results"); mkdir("results"); end
        save_name = sprintf("results/%s", save_name);
    end
    if opt.save_add_time; save_name = sprintf("%s_%s", save_name, string(datetime('now', 'Format', 'y-M-d_HH-mm-ss'))); end

    % Generate test_info structure, which contains info related to the test setup
    test_info.type = 'grid';
    test_info.param = param;
    test_info.seed = seed;
    test_info.n_rows = n_rows;
    test_info.n_cols = n_cols;
    test_info.num_tasks = num_tasks_tt;
    test_info.n = prob.n;
    test_info.n_e = prob.n_e;
    test_info.m = prob.m;
    test_info.k = prob.k;
    test_info.nz = prob.nz;
    test_info.a_loc = a_loc;
    test_info.t_loc = t_loc;
    test_info.solved_heur = false;
    test_info.solved_central = false;

    % Save info structure
    if save_me
        save(save_name, 'test_info');
    end

    %% Solve MRTA problem with PReDuaL

    if run_heur

        fprintf("  > Starting dual Lagrangian heuristic solver at %s.", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')));
        heur_opt = opt.heur;
        heur_opt.solver = opt.solver;
        heur_args = namedargs2cell(heur_opt);
        [~, V_d, info_d] = heur_distLag_MRTA(prob, a_loc, t_loc, heur_args{:});
        fprintf(" Finished at %s with exit flag %d", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')), info_d.e_flag);
        if isfinite(info_d.gap)
            fprintf(" and certified optimality gap %.2f%%", 100*info_d.gap);
        end
        fprintf("\n");
        test_info.solved_heur = true;

        % Generate res structure, containing the results obtained from the solver
        res.time_d_par_total = info_d.run_time_par_total; % Computation time of Alg. 1 + primal recovery (Section VII-B)
        res.time_d_alg1 = info_d.run_time_par_total_alg1; % Computation time of Alg. 1 only
        res.e_flag_alg1_d = info_d.e_flag_alg1; % Exit flag of Alg. 1 (before the primal recovery)
        res.V_alg1_d = info_d.V_alg1; % Cost of the solution of Alg. 1 (before the primal recovery)
        res.LB_d = info_d.LB; % Lower bound of the optimal cost computed by the column generation
        res.gap_d = info_d.gap; % Certified optimality gap (V_d - LB_d)/LB_d
        res.time_d_par = info_d.run_time_par;
        res.time_d_overhead = info_d.run_time_overhead;
        res.time_d_build = info_d.run_time_build;
        res.time_d_all = info_d.run_time_all; % Measured wall-clock time of Alg. 1 (sequential or parallel)
        res.time_d_wall = info_d.run_time_wall_total; % Measured wall-clock time of Alg. 1 + primal recovery
        res.parallel_d = info_d.parallel; % true if the MILP problems of the agents were solved in parallel
        res.e_flag_d = info_d.e_flag;
        res.feasible_d = info_d.feasible;
        res.k_d = info_d.k;
        res.num_inner_tl_d = info_d.num_inner_tl + info_d.num_inner_nosol; % MILPs of Alg. 1 stopped by 'inner_max_time'
        res.V_d = V_d;
        % Computation time of the primal recovery: route selection problems (coordinator), MILP problems
        % of the agents, and the rest (route search, etc.)
        res.time_d_rec = NaN;
        res.time_d_rec_sel = NaN;
        res.time_d_rec_agents = NaN;
        res.time_d_rec_overhead = NaN;
        if isfield(info_d, 'recovery')
            R = info_d.recovery;
            res.time_d_rec = R.run_time_par_total;
            res.time_d_rec_sel = R.run_time_coord;
            res.time_d_rec_agents = R.run_time_par - R.run_time_coord;
            res.time_d_rec_overhead = R.run_time_overhead;
        end

        % Update saved info
        if save_me
            save(save_name, '-append', 'test_info', 'res');
        end

    end

    %% Solve MRTA problem using centralized MILP solver

    if run_central

        fprintf("  > Starting central solver at %s.", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')));
        prob.build_MILP(a_loc, t_loc, 'sparse', true, 'condense', true, 'tight', opt.central_tight); % Build MILP problem (9)
        solver_opt = opt.solver;
        solver_opt.max_time = opt.central_max_time;
        solver_args = namedargs2cell(solver_opt);
        [~, V_c, info_c] = prob.miqp.solve(solver_args{:});
        fprintf(" Finished at %s with exit flag %d", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')), info_c.e_flag);
        test_info.solved_central = info_c.e_flag == 1;

        % Update res structure with results of centralized solver
        res.time_c = info_c.run_time;
        res.e_flag_c = info_c.e_flag;
        res.feasible_c = info_c.feasible;
        res.V_c = V_c;
        res.LB_c = info_c.obj_bound; % Lower bound of the optimal cost given by the solver (at the time limit if not solved)
        res.gap_c = NaN; % Relative gap (V_c - LB_c)/|LB_c| of the central solver (about 0 if solved)
        if isfinite(V_c) && isfinite(res.LB_c)
            res.gap_c = (V_c - res.LB_c)/max(abs(res.LB_c), 1e-12);
        end
        if info_c.e_flag ~= 1 && isfinite(res.gap_c)
            fprintf(" (not solved: remaining gap %.2f%%)", 100*res.gap_c);
        end
        fprintf("\n");

        % Update saved info
        if save_me
            save(save_name, '-append', 'test_info', 'res');
        end

    end

    %% Results if both solvers were run

    if run_central && run_heur

        res.err_V = NaN;
        res.err_V_alg1 = NaN;
        if res.e_flag_c == 1 && res.e_flag_d >= 0
            res.err_V = (res.V_d - res.V_c)/res.V_c; % Cost error (in % in Table I)
            fprintf("  > Time centralized: %g. Time PReDuaL: %g. Cost centralized: %g. Cost PReDuaL: %g. Error: %3.2f%%\n", ...
                res.time_c, res.time_d_par_total, res.V_c, res.V_d, 100*res.err_V);
        end
        if res.e_flag_c == 1 && res.e_flag_alg1_d >= 0
            res.err_V_alg1 = (res.V_alg1_d - res.V_c)/res.V_c; % Cost error of Alg. 1 without primal recovery
            fprintf("  > Without primal recovery: time Alg. 1: %g. Cost Alg. 1: %g. Error: %3.2f%%\n", ...
                res.time_d_alg1, res.V_alg1_d, 100*res.err_V_alg1);
        end
        if save_me
            save(save_name, '-append', 'test_info', 'res');
        end

    end

    fprintf("\n");

end
