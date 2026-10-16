%% analyze_GA_batch() - Compares the GA with PReDuaL and the central MILP on the scenarios of a grid (Table II)
%
% Reads the results of run_GA_batch (results/ga_<name>.mat) and of run_grid_batch (results/<name>.mat) for the
% same scenarios, and prints, for each scenario and on average, the cost errors with respect to the optimal
% cost of the central MILP (9) and the computation times of:
%   - the central MILP (9)
%   - PReDuaL (time of Section VII-B, and measured wall-clock time)
%   - Algorithm 1 alone (before the primal recovery)
%   - the GA (wall-clock time, single thread)
%
% INPUTS:
% @get_n: Number of locations n of the scenarios (e.g., 56 for the 8 x 7 grid). NaN for all scenarios.
%
% NAME/VALUE INPUTS:
% 'folder': Folder with the .mat files (default "results")
% 'time_budget': Results of the GA to compare (see run_GA_batch): "none" (default; GA with its own stopping
%                rule, results/ga_*.mat), "paper" or "wall" (GA with the same time as PReDuaL,
%                results/ga-paper_*.mat or results/ga-wall_*.mat)
%
% OUTPUTS:
% @stats: Structure with the statistics [avrg, median, max, min] (cost errors in %, times in s) and the table
%         of results of each scenario (stats.table)
%
% EXAMPLES: run_grid_batch(8, 7); run_GA_batch(8, 7); analyze_GA_batch(56);
%           run_GA_batch(8, 7, time_budget="paper"); analyze_GA_batch(56, time_budget="paper");
%

function stats = analyze_GA_batch(get_n, opt)
    arguments
        get_n (1, 1) double = NaN
        opt.folder {mustBeText} = "results"
        opt.time_budget {mustBeText} = "none"
    end
    folder = char(opt.folder);
    budget = lower(char(opt.time_budget));
    switch budget
        case 'none'
            pattern = 'ga_*.mat';
            ga_label = 'GA (own stopping rule)';
        case {'paper', 'wall'}
            pattern = ['ga-', budget, '_*.mat'];
            ga_label = sprintf('GA with the time of PReDuaL (%s time)', budget);
        otherwise
            error("analyze_GA_batch:time_budget", 'time_budget must be "none", "paper" or "wall"');
    end

    files = dir(fullfile(folder, pattern));
    names = {};
    % Columns: tasks, V_c, time_c, solved_c, V_d, time_d, wall_d, V_alg1, time_alg1, V_ga, time_ga, V_ga_mean, feasible_ga,
    %          GA runs (within the time budget), V of the GA with its own stopping rule (only with a time budget)
    T = zeros(0, 15);
    for f = 1:length(files)
        data = load(fullfile(folder, files(f).name));
        if ~isfield(data, 'res_ga'); continue; end
        g = data.res_ga;
        if ~isnan(get_n) && g.n ~= get_n; continue; end
        row = nan(1, 15);
        row(1) = g.num_tasks;
        row([10, 11, 12, 13]) = [g.V, g.time, g.V_mean, double(g.feasible)];
        row(14) = get_field(g, 'ga_runs');
        default_file = fullfile(folder, ['ga_', g.name, '.mat']);
        if ~strcmp(budget, 'none') && isfile(default_file)
            d = load(default_file);
            if isfield(d, 'res_ga'); row(15) = d.res_ga.V; end
        end
        res_file = fullfile(folder, [g.name, '.mat']);
        if isfile(res_file)
            saved = load(res_file);
            if isfield(saved, 'res')
                r = saved.res;
                row(2) = get_field(r, 'V_c');
                row(3) = get_field(r, 'time_c');
                row(4) = double(get_field(r, 'e_flag_c') == 1);
                row(5) = get_field(r, 'V_d');
                row(6) = get_field(r, 'time_d_par_total');
                row(7) = get_field(r, 'time_d_wall');
                row(8) = get_field(r, 'V_alg1_d');
                row(9) = get_field(r, 'time_d_alg1');
                if get_field(r, 'e_flag_d') < 0; row(5) = inf; end
                if get_field(r, 'e_flag_alg1_d') < 0; row(8) = inf; end
            end
        end
        names{end+1} = g.name; %#ok<AGROW>
        T(end+1, :) = row; %#ok<AGROW>
    end
    if isempty(names)
        error("analyze_GA_batch:no_tests", "No results of run_GA_batch (%s) in %s/ satisfy the criteria (n = %g).", pattern, folder, get_n);
    end

    solved_c = T(:, 4) == 1 & isfinite(T(:, 2)) & T(:, 2) > 0;
    err = @(V) 100*(V - T(:, 2))./T(:, 2);
    err_d = err(T(:, 5));
    err_alg1 = err(T(:, 8));
    err_ga = err(T(:, 10));
    err_ga_mean = err(T(:, 12));

    %% Table of results of each scenario
    fprintf("*** COMPARISON OF THE GA WITH PReDuaL AND THE CENTRAL MILP (9) ***\n");
    fprintf("*** %s ***\n", ga_label);
    fprintf("%-24s %5s | %10s %8s | %10s %8s %8s | %8s %8s | %10s %8s %8s\n", "Scenario", "Tasks", ...
        "Central", "Time[s]", "PReDuaL", "Err[%]", "Time[s]", "Alg.1[%]", "Time[s]", "GA", "Err[%]", "Time[s]");
    for q = 1:length(names)
        e = [err_d(q), err_alg1(q), err_ga(q)];
        e(~solved_c(q) | isnan(e)) = NaN;
        fprintf("%-24s %5d | %10.4f %8.3f | %10.4f %8.2f %8.3f | %8.2f %8.3f | %10.4f %8.2f %8.3f\n", names{q}, T(q, 1), ...
            T(q, 2), T(q, 3), T(q, 5), e(1), T(q, 6), e(2), T(q, 9), T(q, 10), e(3), T(q, 11));
    end

    %% Statistics
    stats.num_tests = length(names);
    stats.num_central = nnz(solved_c);
    stats.time_c = stat4(T(solved_c, 3));
    stats.time_d = stat4(T(~isnan(T(:, 6)), 6));
    stats.wall_d = stat4(T(~isnan(T(:, 7)), 7));
    stats.time_alg1 = stat4(T(~isnan(T(:, 9)), 9));
    stats.time_ga = stat4(T(:, 11));
    stats.err_d = stat4(err_d(solved_c & isfinite(err_d)));
    stats.err_alg1 = stat4(err_alg1(solved_c & isfinite(err_alg1)));
    stats.err_ga = stat4(err_ga(solved_c & isfinite(err_ga)));
    stats.err_ga_mean = stat4(err_ga_mean(solved_c & isfinite(err_ga_mean)));
    err_ga_default = err(T(:, 15));
    stats.err_ga_default = stat4(err_ga_default(solved_c & isfinite(err_ga_default)));
    stats.feas_ga = nnz(T(:, 13) == 1);
    stats.ga_runs = stat4(T(~isnan(T(:, 14)), 14)); % Runs of the GA (within the time budget)
    has_d = isfinite(T(:, 5)) & isfinite(T(:, 10));
    tol = 1e-6*max(1, abs(T(:, 5))); % Relative tolerance (the MILP solvers stop within a small optimality gap)
    stats.ga_better = nnz(has_d & T(:, 10) < T(:, 5) - tol); % Scenarios where the GA is better than PReDuaL
    stats.ga_equal = nnz(has_d & abs(T(:, 10) - T(:, 5)) <= tol);
    stats.ga_worse = nnz(has_d & T(:, 10) > T(:, 5) + tol);
    stats.names = names;
    stats.table = T;

    fprintf("\n> Scenarios: %d (central MILP solved in %d; the cost errors are computed for these)\n", stats.num_tests, stats.num_central);
    fprintf("                                        %10s %10s %10s %10s\n", "Avrg.", "Med.", "Max.", "Min.");
    fprintf("> Comp. time central MILP (9) [s]       %10.3f %10.3f %10.3f %10.3f\n", stats.time_c);
    fprintf("> Comp. time PReDuaL [s]                %10.3f %10.3f %10.3f %10.3f   (Section VII-B)\n", stats.time_d);
    if ~all(isnan(stats.wall_d))
        fprintf("> Wall-clock PReDuaL [s]                %10.3f %10.3f %10.3f %10.3f\n", stats.wall_d);
    end
    fprintf("> Comp. time Alg. 1 [s]                 %10.3f %10.3f %10.3f %10.3f   (without recovery)\n", stats.time_alg1);
    fprintf("> Comp. time GA [s]                     %10.3f %10.3f %10.3f %10.3f\n", stats.time_ga);
    if ~strcmp(budget, 'none') && ~all(isnan(stats.ga_runs))
        fprintf("> Runs of the GA within the time        %10.1f %10.1f %10d %10d\n", stats.ga_runs);
    end
    fprintf("> Cost error PReDuaL (%%)                %10.2f %10.2f %10.2f %10.2f\n", stats.err_d);
    fprintf("> Cost error Alg. 1 (%%)                 %10.2f %10.2f %10.2f %10.2f   (without recovery)\n", stats.err_alg1);
    fprintf("> Cost error GA (%%)                     %10.2f %10.2f %10.2f %10.2f\n", stats.err_ga);
    if any(abs(T(:, 12) - T(:, 10)) > 1e-12)
        fprintf("> Cost error GA, average of runs (%%)    %10.2f %10.2f %10.2f %10.2f\n", stats.err_ga_mean);
    end
    if ~all(isnan(stats.err_ga_default))
        fprintf("> Cost error GA, own stopping rule (%%)  %10.2f %10.2f %10.2f %10.2f   (results/ga_*.mat)\n", stats.err_ga_default);
    end
    fprintf("> GA: feasible in %d of %d scenarios. Compared with PReDuaL: better in %d, equal in %d, worse in %d\n", ...
        stats.feas_ga, stats.num_tests, stats.ga_better, stats.ga_equal, stats.ga_worse);

end

function val = get_field(s, name)
    % Returns s.(name), or NaN if the field does not exist
    if isfield(s, name) && ~isempty(s.(name))
        val = double(s.(name));
    else
        val = NaN;
    end
end

function s = stat4(x)
    % Returns [average, median, max, min] of x (NaN if x is empty)
    x = x(:);
    if isempty(x)
        s = NaN(1, 4);
    else
        s = [mean(x), median(x), max(x), min(x)];
    end
end
