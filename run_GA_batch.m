%% run_GA_batch() - Solves the random scenarios of run_grid_batch with the GA (Section VII-C)
%
% Generates the same scenarios as run_grid_batch(n_rows, n_cols), solves them with the GA (ga_MRTA),
% and saves the results in results/ga_<save_name>_gs_<n>_ss_<state seed>.mat (structure res_ga).
% Use analyze_GA_batch(n_rows*n_cols) to compare them with the results of run_grid_batch.
%
% With the option 'time_budget', the GA gets the same computation time as PReDuaL on each scenario
% (read from the results of run_grid_batch). The GA restarts with a new initial population each time a
% run stops (50 generations without improvement) until the time is used, and the best solution is
% reported. These results are saved in results/ga-<time_budget>_<save_name>_gs_*.mat.
%
% INPUTS:
% @n_rows: Number of rows of the grid (n_v)
% @n_cols: Number of columns of the grid (m_v)
%
% NAME/VALUE INPUTS:
% 'num_tests': Number of scenarios (default 10, as run_grid_batch)
% 'save_me': Save the results in results/ (default true)
% 'save_name': Prefix of the files of run_grid_batch (default "grid")
% 'runs': Number of runs of the GA on each scenario, with seeds seed, seed + 1, ... (default 1). The best
%         run is reported (res_ga.V, res_ga.time), as well as the averages over the runs (V_mean, time_mean)
% 'ga': Structure with name/value options for ga_MRTA, e.g., struct('pop', 200, 'stop', 200, 'seed', 1)
% 'time_budget': "none" (default): the GA stops with its own rule (200 generations without improvement).
%                "paper": computation time of PReDuaL of Section VII-A (res.time_d_par_total).
%                "wall": measured wall-clock time of PReDuaL (res.time_d_wall).
% 'skip_done': Skip the scenarios whose GA results already exist in results/ (default false). Used to
%              resume an interrupted batch
%
% OUTPUTS:
% @results: Cell array with the structure res_ga of each scenario
%

function results = run_GA_batch(n_rows, n_cols, opt)
    arguments
        n_rows (1, 1) {mustBePositive, mustBeInteger}
        n_cols (1, 1) {mustBePositive, mustBeInteger}
        opt.num_tests (1, 1) {mustBePositive, mustBeInteger} = 10
        opt.save_me (1, 1) {mustBeNumericOrLogical} = true
        opt.save_name {mustBeText} = "grid"
        opt.runs (1, 1) {mustBePositive, mustBeInteger} = 1
        opt.ga struct = struct()
        opt.time_budget {mustBeText} = "none"
        opt.skip_done (1, 1) {mustBeNumericOrLogical} = false
    end
    budget = lower(char(opt.time_budget));
    if ~any(strcmp(budget, {'none', 'paper', 'wall'}))
        error("run_GA_batch:time_budget", 'time_budget must be "none", "paper" or "wall"');
    end
    budget_field = struct('paper', 'time_d_par_total', 'wall', 'time_d_wall');
    if strcmp(budget, 'none')
        prefix = 'ga_';
    else
        prefix = ['ga-', budget, '_'];
    end

    % Same scenarios as run_grid_batch
    n = n_rows * n_cols;
    seed.graph = n;
    seed_state = grid_seeds(seed.graph, opt.num_tests);
    param = gen_param_grid(n_rows, n_cols, n_rows, n_cols);

    ga_opt = opt.ga;
    if ~isfield(ga_opt, 'seed'); ga_opt.seed = 1; end
    first_seed = ga_opt.seed;
    if opt.save_me && ~isfolder('results'); mkdir('results'); end

    results = cell(1, opt.num_tests);
    for s = 1:opt.num_tests
        seed.state = seed_state(s);
        name = sprintf('%s_gs_%d_ss_%d', char(opt.save_name), seed.graph, seed.state);
        if opt.skip_done
            ga_file = fullfile('results', [prefix, name, '.mat']);
            if isfile(ga_file)
                try
                    S = load(ga_file, 'res_ga');
                    results{s} = S.res_ga;
                    fprintf("*** GA: test #%d of %d (%d x %d grid), %s already done: skipped ***\n", s, opt.num_tests, n_rows, n_cols, name);
                    continue
                catch
                    % Corrupted file (e.g., interrupted while saving): solve it again
                end
            end
        end
        [prob, a_loc, t_loc, gen_info] = gen_rand_grid_instance(param, seed);
        fprintf("*** GA: test #%d of %d (%d x %d grid), %s: m = %d robots, %d tasks ***\n", s, opt.num_tests, ...
            n_rows, n_cols, name, prob.m, gen_info.num_tasks);

        % Results of run_grid_batch for the same scenario (if available)
        saved = struct();
        res_file = fullfile('results', [name, '.mat']);
        if isfile(res_file)
            saved = load(res_file);
            if isfield(saved, 'test_info') && (~isequal(saved.test_info.a_loc, a_loc) || ~isequal(saved.test_info.t_loc, t_loc))
                warning("run_GA_batch:scenario", "The scenario of %s differs from the one generated now", res_file);
                saved = struct();
            end
        end

        % Time budget of the GA (time of PReDuaL on this scenario)
        T_budget = NaN;
        ga_opt_s = ga_opt;
        if ~strcmp(budget, 'none')
            if isfield(saved, 'res') && isfield(saved.res, budget_field.(budget))
                T_budget = double(saved.res.(budget_field.(budget)));
            end
            if ~(isfinite(T_budget) && T_budget > 0)
                warning("run_GA_batch:time_budget", "Scenario %s skipped: no time of PReDuaL in %s (run run_grid_batch first)", name, res_file);
                fprintf("\n");
                continue
            end
            ga_opt_s.time = T_budget;
            if ~isfield(ga_opt_s, 'restarts'); ga_opt_s.restarts = true; end
            % With a time budget, restart after 50 generations without improvement (unless 'stop' is given)
            if ~isfield(ga_opt_s, 'stop'); ga_opt_s.stop = 50; end
        end

        % Runs of the GA
        V_runs = inf(1, opt.runs);
        time_runs = nan(1, opt.runs);
        time_best_runs = nan(1, opt.runs);
        gens_runs = nan(1, opt.runs);
        seeds = zeros(1, opt.runs);
        best = [];
        for r = 1:opt.runs
            ga_opt_s.seed = first_seed;
            if first_seed > 0; ga_opt_s.seed = first_seed + r - 1; end
            ga_opt_s.name = name;
            seeds(r) = ga_opt_s.seed;
            ga_args = namedargs2cell(ga_opt_s);
            [~, V, info] = ga_MRTA(prob, a_loc, t_loc, ga_args{:});
            V_runs(r) = V;
            time_runs(r) = info.time;
            time_best_runs(r) = info.time_best;
            gens_runs(r) = info.generations;
            if isempty(best) || V < best.V
                best = struct('V', V, 'info', info);
            end
        end

        res_ga = struct();
        res_ga.name = name;
        res_ga.seed = seed;
        res_ga.n_rows = gen_info.n_rows;
        res_ga.n_cols = gen_info.n_cols;
        res_ga.n = prob.n;
        res_ga.m = prob.m;
        res_ga.k = prob.k;
        res_ga.num_tasks = gen_info.num_tasks;
        res_ga.V = best.V; % Cost (8) of the best solution of the GA (Inf if no feasible solution was found)
        res_ga.e_flag = best.info.e_flag; % 1: feasible solution of (9); -1: no feasible solution
        res_ga.feasible = best.info.feasible;
        res_ga.time = best.info.time; % Computation time of the GA [s] (best run)
        res_ga.time_best = best.info.time_best; % Time at which the GA found its best solution [s] (best run)
        res_ga.wall_time = best.info.wall_time; % Wall-clock time of the call to the GA executable [s] (best run)
        res_ga.generations = best.info.generations;
        res_ga.ga_runs = best.info.runs; % Runs of the GA within the time budget (1 without time budget)
        res_ga.budget = budget; % 'none', 'paper' or 'wall'
        res_ga.time_budget = T_budget; % Time budget of the GA [s] (NaN if none)
        res_ga.V_mean = mean(V_runs); % Averages over the runs
        res_ga.time_mean = mean(time_runs);
        res_ga.runs = opt.runs;
        res_ga.seeds = seeds;
        res_ga.V_runs = V_runs;
        res_ga.time_runs = time_runs;
        res_ga.time_best_runs = time_best_runs;
        res_ga.generations_runs = gens_runs;
        res_ga.routes = best.info.routes; % Route (locations) and tasks of each robot (best run)
        res_ga.tasks = best.info.tasks;
        res_ga.task_v = best.info.task_v; % Location and type of each task
        res_ga.task_t = best.info.task_t;
        res_ga.ga_opt = rmfield(ga_opt_s, {'seed', 'name'});

        % Print the results (and the comparison with the results of run_grid_batch)
        line = sprintf('  > GA: cost %.6g, time %.3f s (best solution at %.3f s), %d generations', res_ga.V, res_ga.time, ...
            res_ga.time_best, res_ga.generations);
        if ~strcmp(budget, 'none')
            line = [line, sprintf(' in %d runs (time budget %.3f s, %s)', res_ga.ga_runs, T_budget, budget)];
        end
        if opt.runs > 1
            line = [line, sprintf(' (best of %d runs; average cost %.6g, time %.3f s)', opt.runs, res_ga.V_mean, res_ga.time_mean)];
        end
        if ~res_ga.feasible
            line = [line, sprintf(' (no feasible solution: %s)', best.info.msg)];
        end
        fprintf("%s\n", line);
        if isfield(saved, 'res') && isfield(saved.res, 'V_c') && saved.res.e_flag_c == 1
            line = sprintf('  > Central: cost %.6g, time %.3f s. GA error: %.2f%%', saved.res.V_c, saved.res.time_c, ...
                100*(res_ga.V - saved.res.V_c)/saved.res.V_c);
            if isfield(saved.res, 'V_d')
                line = [line, sprintf('. PReDuaL error: %.2f%%', 100*(saved.res.V_d - saved.res.V_c)/saved.res.V_c)];
            end
            fprintf("%s\n", line);
        end
        fprintf("\n");

        if opt.save_me
            save(fullfile('results', [prefix, name, '.mat']), 'res_ga');
        end
        results{s} = res_ga;
    end

end
