%% analyze_rand_batch() - Summarizes the results of run_grid_batch() as in Table I of the article
%
% INPUTS:
% @get_n: Number of locations n of the tests to analyze (e.g., 56 for the 8 x 7 grid). NaN for all tests.
% @folder: Folder containing the .mat files (default "results")
%
% OUTPUTS:
% @stats: Structure containing the statistics [avrg, median, max, min] of the results
%
% EXAMPLE: analyze_rand_batch(56)
%

function stats = analyze_rand_batch(get_n, folder)
    arguments
        get_n (1, 1) double = NaN
        folder {mustBeText} = "results"
    end

    %% Load the results of eligible tests
    mat = dir(fullfile(folder, '*.mat'));
    data = struct('nz', {}, 'grid', {}, 'time_c', {}, 'e_flag_c', {}, 'time_d', {}, 'e_flag_d', {}, 'k_d', {}, 'err_V', {}, ...
                  'time_alg1', {}, 'e_flag_alg1', {}, 'err_V_alg1', {}, 'gap', {}, 'wall_d', {}, 'wall_alg1', {}, 'parallel', {}, ...
                  'rec', {}, 'rec_sel', {}, 'rec_agents', {}, 'rec_ovh', {}, 'V_c', {}, 'LB_c', {}, 'gap_c', {}, 'V_d', {}, 'LB_d', {});
    for i = 1:length(mat)
        test = load(fullfile(mat(i).folder, mat(i).name));
        if ~isfield(test, 'test_info') || ~isfield(test, 'res') || ~isfield(test.test_info, 'n_rows')
            continue; % Not a result of run_rand_grid
        end
        if ~isnan(get_n) && test.test_info.n ~= get_n
            continue;
        end
        d.nz = test.test_info.nz;
        d.grid = sprintf('%d x %d', test.test_info.n_rows, test.test_info.n_cols);
        d.time_c = get_field(test.res, 'time_c');
        d.e_flag_c = get_field(test.res, 'e_flag_c');
        d.time_d = get_field(test.res, 'time_d_par_total');
        d.e_flag_d = get_field(test.res, 'e_flag_d');
        d.k_d = get_field(test.res, 'k_d');
        d.err_V = get_field(test.res, 'err_V');
        d.time_alg1 = get_field(test.res, 'time_d_alg1');
        d.e_flag_alg1 = get_field(test.res, 'e_flag_alg1_d');
        d.err_V_alg1 = get_field(test.res, 'err_V_alg1');
        d.gap = get_field(test.res, 'gap_d');
        d.wall_d = get_field(test.res, 'time_d_wall');
        d.wall_alg1 = get_field(test.res, 'time_d_all');
        d.parallel = get_field(test.res, 'parallel_d');
        d.rec = get_field(test.res, 'time_d_rec'); % Breakdown of the time of the primal recovery
        d.rec_sel = get_field(test.res, 'time_d_rec_sel');
        d.rec_agents = get_field(test.res, 'time_d_rec_agents');
        d.rec_ovh = get_field(test.res, 'time_d_rec_overhead');
        d.V_c = get_field(test.res, 'V_c');
        d.LB_c = get_field(test.res, 'LB_c'); % Lower bound of the central solver
        d.gap_c = get_field(test.res, 'gap_c');
        d.V_d = get_field(test.res, 'V_d');
        d.LB_d = get_field(test.res, 'LB_d');
        data(end+1) = d; %#ok<AGROW>
    end

    num_tests = length(data);
    if num_tests == 0
        error("analyze_rand_batch:no_tests", "No tests in %s/ satisfy the criteria (n = %g).", folder, get_n);
    end

    %% Compute statistics
    e_flag_c = [data.e_flag_c];
    e_flag_d = [data.e_flag_d];
    solved_c = e_flag_c == 1;
    solved_d = ~isnan(e_flag_d);
    time_c = [data.time_c];
    time_d = [data.time_d];
    k_d = [data.k_d];
    err_V = 100*[data.err_V];
    time_alg1 = [data.time_alg1];
    wall_d = [data.wall_d];
    wall_alg1 = [data.wall_alg1];
    e_flag_alg1 = [data.e_flag_alg1];
    err_V_alg1 = 100*[data.err_V_alg1];
    gap = 100*[data.gap];

    stats.num_tests = num_tests;
    stats.grid = unique({data.grid});
    stats.nz = round(mean([data.nz]));
    stats.time_c = stat4(time_c(solved_c)); % Comp. time (10) [s]
    stats.time_d = stat4(time_d(solved_d)); % Comp. time PReDuaL [s] (Section VII-B)
    stats.wall_d = stat4(wall_d(solved_d & ~isnan(wall_d))); % Measured wall-clock time PReDuaL [s]
    stats.wall_alg1 = stat4(wall_alg1(solved_d & ~isnan(wall_alg1))); % Measured wall-clock time Alg. 1 [s]
    stats.parallel = nnz([data.parallel] == 1); % Number of tests solved with parallel MILP solves
    stats.k_d = stat4(k_d(solved_d)); % Iterations Alg. 1
    stats.err_V = stat4(err_V(~isnan(err_V))); % Cost error (%)
    stats.solved_c = nnz(solved_c); % Number of central MILPs (10) solved
    stats.feas_d = nnz(e_flag_d == 1); % Number of feasible solutions of Alg. 1
    stats.over_d = nnz(e_flag_d == 2 | e_flag_d == 0); % Number of over-feasible solutions of Alg. 1
    stats.infeas_d = nnz(e_flag_d == -1); % Number of tests where Alg. 1 did not find an (over-)feasible solution
    stats.gap = stat4(gap(~isnan(gap))); % Optimality gap (%) from the column generation
    stats.certified = nnz(gap <= 0.1); % Number of solutions certified to be within 0.1% of the optimal cost
    % Results of Algorithm 1 alone (before the primal recovery)
    stats.time_alg1 = stat4(time_alg1(solved_d & ~isnan(time_alg1)));
    stats.err_V_alg1 = stat4(err_V_alg1(~isnan(err_V_alg1)));
    stats.feas_alg1 = nnz(e_flag_alg1 == 1);
    stats.over_alg1 = nnz(e_flag_alg1 == 2 | e_flag_alg1 == 0);
    % Tests where the central MILP was not solved (e.g., time limit): remaining gap of the central solver,
    % comparison of PReDuaL with its best solution, and gap of PReDuaL with respect to the best lower bound
    unsolved_c = ~solved_c & solved_d;
    V_c = [data.V_c];
    V_d = [data.V_d];
    LB_best = max([data.LB_c], [data.LB_d]); % (max ignores NaN)
    gap_c = 100*[data.gap_c];
    gap_best = 100*(V_d - LB_best)./abs(LB_best);
    tol = 1e-6*max(1, abs(V_c));
    stats.unsolved_c = nnz(unsolved_c);
    stats.gap_c = stat4(gap_c(unsolved_c & isfinite(gap_c))); % Remaining gap of the central solver (%)
    stats.d_vs_c = [nnz(unsolved_c & V_d < V_c - tol), nnz(unsolved_c & abs(V_d - V_c) <= tol), ...
                    nnz(unsolved_c & V_d > V_c + tol), nnz(unsolved_c & ~isfinite(V_c))]; % better / equal / worse / no central solution
    stats.gap_best = stat4(gap_best(unsolved_c & isfinite(gap_best))); % Gap of PReDuaL w.r.t. the best bound (%)

    %% Print results
    fprintf("*** REPORT (Table I) ***\n");
    fprintf("> Number of loaded tests: %d (n: %g, grid: %s)\n", num_tests, get_n, strjoin(stats.grid, ', '));
    fprintf("> Decision variables nz (avrg.): %d\n", stats.nz);
    fprintf("                                   %10s %10s %10s %10s\n", "Avrg.", "Med.", "Max.", "Min.");
    fprintf("> Comp. time (10) [s]              %10.2f %10.2f %10.2f %10.2f\n", stats.time_c);
    fprintf("> Comp. time PReDuaL [s]           %10.2f %10.2f %10.2f %10.2f\n", stats.time_d);
    if ~all(isnan(stats.wall_d))
        fprintf("> Wall-clock PReDuaL [s]           %10.2f %10.2f %10.2f %10.2f   (parallel in %d of %d tests)\n", ...
            stats.wall_d, stats.parallel, num_tests);
    end
    rec = [data.rec];
    has_rec = solved_d & ~isnan(rec) & ~isnan(time_alg1);
    if any(has_rec)
        % Average breakdown of the computation time of PReDuaL
        stats.time_breakdown = [mean(time_alg1(has_rec)), mean(rec(has_rec)), mean([data(has_rec).rec_sel]), ...
                                mean([data(has_rec).rec_agents]), mean([data(has_rec).rec_ovh])];
        fprintf("  > Avrg. time: Alg. 1 %.2f s, recovery %.2f s (route selection %.2f s, agents %.2f s, route search and other %.2f s)\n", ...
            stats.time_breakdown);
    end
    fprintf("> Iterations Alg. 1                %10.1f %10.1f %10d %10d\n", stats.k_d);
    fprintf("> Cost error (%%)                   %10.2f %10.2f %10.2f %10.2f\n", stats.err_V);
    if ~all(isnan(stats.gap))
        fprintf("> Optimality gap (%%)               %10.2f %10.2f %10.2f %10.2f   (within 0.1%%: %d)\n", stats.gap, stats.certified);
    end
    fprintf("> Sol. (10) #solved: %d\n", stats.solved_c);
    if stats.unsolved_c > 0
        fprintf("> Central MILP (10) not solved in %d tests (e.g., time limit):\n", stats.unsolved_c);
        if ~all(isnan(stats.gap_c))
            fprintf("  > Remaining gap of the central solver (%%)        %10.2f %10.2f %10.2f %10.2f\n", stats.gap_c);
        end
        fprintf("  > PReDuaL vs. best central solution: better %d, equal %d, worse %d, no central solution %d\n", stats.d_vs_c);
        if ~all(isnan(stats.gap_best))
            fprintf("  > Gap of PReDuaL w.r.t. best bound (%%)           %10.2f %10.2f %10.2f %10.2f\n", stats.gap_best);
        end
    end
    fprintf("> Sol. Alg. 1 (feas. / over): (%d / %d)", stats.feas_d, stats.over_d);
    if stats.infeas_d > 0
        fprintf(", not found: %d", stats.infeas_d);
    end
    fprintf("\n");
    if ~all(isnan(stats.err_V_alg1)) || ~all(isnan(stats.time_alg1))
        fprintf("> Algorithm 1 alone (before the primal recovery):\n");
        fprintf("  > Comp. time Alg. 1 [s] %10.2f %10.2f %10.2f %10.2f\n", stats.time_alg1);
        if ~all(isnan(stats.wall_alg1))
            fprintf("  > Wall-clock Alg. 1 [s] %10.2f %10.2f %10.2f %10.2f\n", stats.wall_alg1);
        end
        fprintf("  > Cost error (%%)        %10.2f %10.2f %10.2f %10.2f\n", stats.err_V_alg1);
        fprintf("  > Sol. Alg. 1 (feas. / over): (%d / %d)\n", stats.feas_alg1, stats.over_alg1);
    end

end

function val = get_field(s, name)
    % Returns s.(name), or NaN if the field does not exist
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = NaN;
    end
end

function s = stat4(x)
    % Returns [average, median, max, min] of x (NaN if x is empty)
    if isempty(x)
        s = NaN(1, 4);
    else
        s = [mean(x), median(x), max(x), min(x)];
    end
end
