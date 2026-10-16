%% run_grid_batch() - Solves the random scenarios of a grid of Table I of the article
%
% Generates an n_rows x n_cols grid with m = n_rows agents and k = n_cols task-types
% (Section VII-B), and solves num_tests scenarios (random agent locations and pending tasks)
% with PReDuaL and the central MILP (9). Results are saved in results/ and can be
% summarized with analyze_rand_batch(n_rows*n_cols).
%
% INPUTS:
% @n_rows: Number of rows of the grid (n_v)
% @n_cols: Number of columns of the grid (m_v)
%
% NAME/VALUE INPUTS:
% 'num_tests': Number of scenarios (default 10)
% 'save_me': Save the results in results/ (default true)
% 'save_name': Prefix of the saved files (default "grid")
% 'ignore_seed_central', 'ignore_seed_heur', 'heur', 'solver', 'central_max_time', 'central_tight': see run_rand_grid
% 'skip_done': Skip the scenarios whose results in results/ are complete, i.e., contain the results of
%              PReDuaL and of the central MILP (default false). Used to resume an interrupted batch
%

function run_grid_batch(n_rows, n_cols, opt)
    arguments
        n_rows (1, 1) {mustBePositive, mustBeInteger}
        n_cols (1, 1) {mustBePositive, mustBeInteger}
        opt.num_tests (1, 1) {mustBePositive, mustBeInteger} = 10
        opt.save_me (1, 1) {mustBeNumericOrLogical} = true
        opt.save_name {mustBeText} = "grid"
        opt.ignore_seed_central = []
        opt.ignore_seed_heur = []
        opt.heur struct = struct()
        opt.solver struct = struct()
        opt.central_max_time (1, 1) double {mustBePositive} = inf
        opt.central_tight (1, 1) {mustBeNumericOrLogical} = false
        opt.skip_done (1, 1) {mustBeNumericOrLogical} = false
    end

    n = n_rows * n_cols; % Number of locations
    m = n_rows; % Number of agents
    k = n_cols; % Number of task-types

    % Random seeds for the agent locations and pending tasks (state) of each scenario
    seed.graph = n; % Use dimension of grid to determine RNG seed
    seed_state = grid_seeds(seed.graph, opt.num_tests);

    % Parameters of the grid
    param = gen_param_grid(n_rows, n_cols, m, k);

    fprintf("*** Using MILP solver: %s ***\n", get_solver_name(opt.solver));

    % Start the parallel pool before the tests (if the MILP problems of the robots are solved in parallel)
    if isfield(opt.heur, 'parallel') && opt.heur.parallel && exist('gcp', 'file') > 0 ...
            && license('test', 'Distrib_Computing_Toolbox')
        try
            pool = gcp();
            fprintf("*** MILP problems of the robots solved in parallel (%d workers) ***\n", pool.NumWorkers);
        catch ME
            warning("run_grid_batch:parallel", "Could not start the parallel pool (%s)", ME.message);
        end
    end

    % Run tests
    for i = 1:opt.num_tests
        seed.state = seed_state(i);
        if opt.skip_done && is_done(sprintf('%s_gs_%d_ss_%d', char(opt.save_name), seed.graph, seed.state))
            fprintf("*** Test #%d of %d (%d x %d grid) already done: skipped ***\n", i, opt.num_tests, n_rows, n_cols);
            continue
        end
        fprintf("*** Running test #%d of %d (%d x %d grid) ***\n", i, opt.num_tests, n_rows, n_cols);
        run_rand_grid(param, seed, opt.save_me, opt.save_name, 'save_add_time', false, ...
            'ignore_seed_central', opt.ignore_seed_central, 'ignore_seed_heur', opt.ignore_seed_heur, ...
            'heur', opt.heur, 'solver', opt.solver, 'central_max_time', opt.central_max_time, ...
            'central_tight', opt.central_tight);
    end

end

function done = is_done(name)
    % true if results/<name>.mat contains the results of PReDuaL and of the central MILP
    done = false;
    file = fullfile('results', [name, '.mat']);
    if ~isfile(file); return; end
    try
        S = load(file, 'res');
        done = isfield(S, 'res') && isfield(S.res, 'V_d') && isfield(S.res, 'e_flag_c');
    catch
        done = false; % Corrupted file (e.g., interrupted while saving)
    end
end

function name = get_solver_name(solver_opt)
    if isfield(solver_opt, 'solver') && solver_opt.solver ~= "auto"
        name = string(solver_opt.solver);
    else
        name = MIQPholder.default_solver();
    end
end
