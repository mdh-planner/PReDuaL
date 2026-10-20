%% run_CBBA_batch() - Solves the random scenarios of run_grid_batch with CBBA
%
% Generates the same scenarios as run_grid_batch(n_rows, n_cols), solves them with the Consensus-Based Bundle
% Algorithm (cbba_MRTA), and saves the results in results/cbba_<save_name>_gs_<n>_ss_<state seed>.mat
% (structure res_cbba). Use analyze_CBBA_batch(n_rows*n_cols) to compare them with PReDuaL, the central MILP
% and the GA.
%
% INPUTS:
% @n_rows: Number of rows of the grid (n_v)
% @n_cols: Number of columns of the grid (m_v)
%
% NAME/VALUE INPUTS:
% 'num_tests': Number of scenarios (default 10, as run_grid_batch)
% 'save_me': Save the results in results/ (default true)
% 'save_name': Prefix of the files of run_grid_batch (default "grid")
% 'cbba': Structure with name/value options for cbba_MRTA, e.g., struct('network', "ring", 'reorder', true)
% 'tag': Name of the variant of CBBA (text without spaces). The results are saved in results/cbba-<tag>_*.mat
%        instead of results/cbba_*.mat, so that several variants can be compared (default "": none)
% 'skip_done': Skip the scenarios whose results already exist in results/ (default false)
%
% OUTPUTS:
% @results: Cell array with the structure res_cbba of each scenario
%

function results = run_CBBA_batch(n_rows, n_cols, opt)
    arguments
        n_rows (1, 1) {mustBePositive, mustBeInteger}
        n_cols (1, 1) {mustBePositive, mustBeInteger}
        opt.num_tests (1, 1) {mustBePositive, mustBeInteger} = 10
        opt.save_me (1, 1) {mustBeNumericOrLogical} = true
        opt.save_name {mustBeText} = "grid"
        opt.cbba struct = struct()
        opt.tag {mustBeText} = ""
        opt.skip_done (1, 1) {mustBeNumericOrLogical} = false
    end
    prefix = cbba_prefix(opt.tag);

    % Same scenarios as run_grid_batch
    n = n_rows * n_cols;
    seed.graph = n;
    seed_state = grid_seeds(seed.graph, opt.num_tests);
    param = gen_param_grid(n_rows, n_cols, n_rows, n_cols);
    if opt.save_me && ~isfolder('results'); mkdir('results'); end

    results = cell(1, opt.num_tests);
    for s = 1:opt.num_tests
        seed.state = seed_state(s);
        name = sprintf('%s_gs_%d_ss_%d', char(opt.save_name), seed.graph, seed.state);
        out_file = fullfile('results', [prefix, name, '.mat']);
        if opt.skip_done && isfile(out_file)
            try
                S = load(out_file, 'res_cbba');
                results{s} = S.res_cbba;
                fprintf("*** CBBA: test #%d of %d (%d x %d grid), %s already done: skipped ***\n", s, opt.num_tests, n_rows, n_cols, name);
                continue
            catch
                % Corrupted file (e.g., interrupted while saving): solve it again
            end
        end
        [prob, a_loc, t_loc, gen_info] = gen_rand_grid_instance(param, seed);
        fprintf("*** CBBA: test #%d of %d (%d x %d grid), %s: m = %d robots, %d tasks ***\n", s, opt.num_tests, ...
            n_rows, n_cols, name, prob.m, gen_info.num_tasks);

        % Results of run_grid_batch for the same scenario (if available)
        saved = struct();
        res_file = fullfile('results', [name, '.mat']);
        if isfile(res_file)
            saved = load(res_file);
            if isfield(saved, 'test_info') && (~isequal(saved.test_info.a_loc, a_loc) || ~isequal(saved.test_info.t_loc, t_loc))
                warning("run_CBBA_batch:scenario", "The scenario of %s differs from the one generated now", res_file);
                saved = struct();
            end
        end

        cbba_opt = opt.cbba;
        cbba_opt.name = name;
        args = namedargs2cell(cbba_opt);
        [~, V, info] = cbba_MRTA(prob, a_loc, t_loc, args{:});

        res_cbba = struct();
        res_cbba.name = name;
        res_cbba.seed = seed;
        res_cbba.n_rows = gen_info.n_rows;
        res_cbba.n_cols = gen_info.n_cols;
        res_cbba.n = prob.n;
        res_cbba.m = prob.m;
        res_cbba.k = prob.k;
        res_cbba.num_tasks = gen_info.num_tasks;
        res_cbba.V = V; % Cost (8) of the solution of CBBA (Inf if it is not feasible)
        res_cbba.e_flag = info.e_flag; % 1: feasible solution of (9); -1: no feasible solution
        res_cbba.feasible = info.feasible;
        res_cbba.converged = info.converged;
        res_cbba.time = info.time; % Serial computation time [s]
        res_cbba.time_par = info.time_par; % Distributed computation time [s] (as the time of PReDuaL)
        res_cbba.time_best = info.time_best; % Serial time at which the final solution was known [s]
        res_cbba.time_par_best = info.time_par_best; % Distributed time at which the final solution was known [s]
        res_cbba.rounds_best = info.rounds_best; % Last round that changed the bundles or winners
        res_cbba.wall_time = info.wall_time; % Wall-clock time of the call to the executable [s]
        res_cbba.rounds = info.rounds;
        res_cbba.messages = info.messages;
        res_cbba.network = info.network;
        res_cbba.diameter = info.diameter;
        res_cbba.routes = info.routes; % Route (locations) and tasks of each robot
        res_cbba.tasks = info.tasks;
        res_cbba.task_v = info.task_v; % Location and type of each task
        res_cbba.task_t = info.task_t;
        res_cbba.tag = char(opt.tag);
        res_cbba.cbba_opt = rmfield(cbba_opt, 'name');

        % Print the results (and the comparison with the results of run_grid_batch)
        line = sprintf('  > CBBA: cost %.6g, %d rounds, time %.4f s distributed (%.4f s serial), final solution after %d rounds (%.4f s distributed)', ...
            res_cbba.V, res_cbba.rounds, res_cbba.time_par, res_cbba.time, res_cbba.rounds_best, res_cbba.time_par_best);
        if ~res_cbba.converged; line = [line, ' (not converged)']; end
        if ~res_cbba.feasible; line = [line, sprintf(' (no feasible solution: %s)', info.msg)]; end
        fprintf("%s\n", line);
        if isfield(saved, 'res') && isfield(saved.res, 'V_c') && saved.res.e_flag_c == 1
            line = sprintf('  > Central: cost %.6g, time %.3f s. CBBA error: %.2f%%', saved.res.V_c, saved.res.time_c, ...
                100*(res_cbba.V - saved.res.V_c)/saved.res.V_c);
            if isfield(saved.res, 'V_d')
                line = [line, sprintf('. PReDuaL error: %.2f%%', 100*(saved.res.V_d - saved.res.V_c)/saved.res.V_c)];
            end
            fprintf("%s\n", line);
        end
        fprintf("\n");

        if opt.save_me
            save(out_file, 'res_cbba');
        end
        results{s} = res_cbba;
    end

end

function prefix = cbba_prefix(tag)
    % Prefix of the result files of a variant of CBBA (see analyze_CBBA_batch)
    tag = char(tag);
    if isempty(tag)
        prefix = 'cbba_';
    elseif any(isspace(tag)) || any(tag == '_')
        error("run_CBBA_batch:tag", "The tag cannot contain spaces or underscores");
    else
        prefix = ['cbba-', tag, '_'];
    end
end
