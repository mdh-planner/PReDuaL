%% analyze_CBBA_batch() - Compares CBBA with PReDuaL, the GA and the central MILP on the scenarios of a grid
%
% Reads the results of run_CBBA_batch (results/cbba_<name>.mat), of run_grid_batch (results/<name>.mat) and,
% if they exist, of run_GA_batch (results/ga-paper_<name>.mat by default) for the same scenarios, and prints,
% for each scenario and on average, the cost errors with respect to the optimal cost of the central MILP (9)
% and the computation times of:
%   - the central MILP (9)
%   - PReDuaL (time of Section VII-A: the MILP problems of the robots are solved in parallel)
%   - Algorithm 1 alone (before the primal recovery)
%   - the GA (single thread)
%   - CBBA (distributed time: the robots compute in parallel, as for PReDuaL; and serial time)
% and the time at which each method found the solution it returns (time to the best solution; NaN for results
% saved before it was recorded).
%
% INPUTS:
% @get_n: Number of locations n of the scenarios (e.g., 56 for the 8 x 7 grid). NaN for all scenarios.
%
% NAME/VALUE INPUTS:
% 'folder': Folder with the .mat files (default "results")
% 'tag': Variant of CBBA (see run_CBBA_batch; default "": results/cbba_*.mat)
% 'ga_budget': Results of the GA to show (see run_GA_batch): "paper" (default, results/ga-paper_*.mat),
%              "wall", "none" (results/ga_*.mat) or "" (do not show the GA)
%
% OUTPUTS:
% @stats: Structure with the statistics [avrg, median, max, min] (cost errors in %, times in s) and the table
%         of results of each scenario (stats.table)
%
% EXAMPLE: run_grid_batch(8, 7); run_GA_batch(8, 7, time_budget="paper"); run_CBBA_batch(8, 7);
%          analyze_CBBA_batch(56);
%

function stats = analyze_CBBA_batch(get_n, opt)
    arguments
        get_n (1, 1) double = NaN
        opt.folder {mustBeText} = "results"
        opt.tag {mustBeText} = ""
        opt.ga_budget {mustBeText} = "paper"
    end
    folder = char(opt.folder);
    if isempty(char(opt.tag))
        prefix = 'cbba_';
    else
        prefix = ['cbba-', char(opt.tag), '_'];
    end
    switch lower(char(opt.ga_budget))
        case ''
            ga_prefix = '';
        case 'none'
            ga_prefix = 'ga_';
        case {'paper', 'wall'}
            ga_prefix = ['ga-', lower(char(opt.ga_budget)), '_'];
        otherwise
            error("analyze_CBBA_batch:ga_budget", 'ga_budget must be "paper", "wall", "none" or ""');
    end

    files = dir(fullfile(folder, [prefix, '*.mat']));
    names = {};
    % Columns: 1 tasks, 2 V_c, 3 time_c, 4 solved_c, 5 V_d, 6 time_d, 7 V_alg1, 8 time_alg1, 9 V_ga, 10 time_ga,
    %          11 V_cbba, 12 time_par_cbba, 13 time_cbba, 14 rounds, 15 converged, 16 feasible_cbba, 17 messages,
    %          time to the best solution: 18 central, 19 PReDuaL, 20 GA, 21 CBBA (distributed), 22 CBBA (serial)
    T = zeros(0, 22);
    for f = 1:length(files)
        data = load(fullfile(folder, files(f).name));
        if ~isfield(data, 'res_cbba'); continue; end
        c = data.res_cbba;
        if ~isnan(get_n) && c.n ~= get_n; continue; end
        row = nan(1, 22);
        row(1) = c.num_tasks;
        row(11:17) = [c.V, c.time_par, c.time, c.rounds, double(c.converged), double(c.feasible), c.messages];
        row(21:22) = [get_field(c, 'time_par_best'), get_field(c, 'time_best')];
        res_file = fullfile(folder, [c.name, '.mat']);
        if isfile(res_file)
            saved = load(res_file);
            if isfield(saved, 'res')
                r = saved.res;
                row(2:8) = [get_field(r, 'V_c'), get_field(r, 'time_c'), double(get_field(r, 'e_flag_c') == 1), ...
                            get_field(r, 'V_d'), get_field(r, 'time_d_par_total'), get_field(r, 'V_alg1_d'), ...
                            get_field(r, 'time_d_alg1')];
                row(18:19) = [get_field(r, 'time_best_c'), get_field(r, 'time_best_d')];
                if get_field(r, 'e_flag_d') < 0; row(5) = inf; end
                if get_field(r, 'e_flag_alg1_d') < 0; row(7) = inf; end
            end
        end
        ga_file = fullfile(folder, [ga_prefix, c.name, '.mat']);
        if ~isempty(ga_prefix) && isfile(ga_file)
            g = load(ga_file);
            if isfield(g, 'res_ga'); row([9, 10, 20]) = [g.res_ga.V, g.res_ga.time, get_field(g.res_ga, 'time_best')]; end
        end
        names{end+1} = c.name; %#ok<AGROW>
        T(end+1, :) = row; %#ok<AGROW>
    end
    if isempty(names)
        error("analyze_CBBA_batch:no_tests", "No results of run_CBBA_batch (%s*.mat) in %s/ satisfy the criteria (n = %g).", prefix, folder, get_n);
    end
    has_ga = any(~isnan(T(:, 9)));

    solved_c = T(:, 4) == 1 & isfinite(T(:, 2)) & T(:, 2) > 0;
    err = @(V) 100*(V - T(:, 2))./T(:, 2);
    err_d = err(T(:, 5));
    err_alg1 = err(T(:, 7));
    err_ga = err(T(:, 9));
    err_cbba = err(T(:, 11));

    %% Table of results of each scenario
    fprintf("*** COMPARISON OF CBBA WITH PReDuaL%s AND THE CENTRAL MILP (9) ***\n", repmat(', THE GA', 1, has_ga));
    fprintf("%-24s %5s | %10s %8s | %8s %8s | %8s | %8s %8s | %10s %8s %8s %6s\n", "Scenario", "Tasks", "Central", ...
        "Time[s]", "PReDuaL%", "Time[s]", "Alg.1%", "GA[%]", "Time[s]", "CBBA", "Err[%]", "Time[s]", "Rounds");
    for q = 1:length(names)
        e = [err_d(q), err_alg1(q), err_ga(q), err_cbba(q)];
        e(~solved_c(q) | isnan(e)) = NaN;
        fprintf("%-24s %5d | %10.4f %8.3f | %8.2f %8.3f | %8.2f | %8.2f %8.3f | %10.4f %8.2f %8.4f %6d%s\n", names{q}, ...
            T(q, 1), T(q, 2), T(q, 3), e(1), T(q, 6), e(2), e(3), T(q, 10), T(q, 11), e(4), T(q, 12), T(q, 14), ...
            repmat(' (not converged)', 1, T(q, 15) == 0));
    end

    %% Statistics
    ok = @(e) e(solved_c & isfinite(e));
    stats.num_tests = length(names);
    stats.num_central = nnz(solved_c);
    stats.time_c = stat4(T(solved_c, 3));
    stats.time_d = stat4(T(~isnan(T(:, 6)), 6));
    stats.time_alg1 = stat4(T(~isnan(T(:, 8)), 8));
    stats.time_ga = stat4(T(~isnan(T(:, 10)), 10));
    stats.time_par_cbba = stat4(T(:, 12));
    stats.time_cbba = stat4(T(:, 13));
    stats.rounds_cbba = stat4(T(:, 14));
    best = @(col, rows) stat4(T(rows & ~isnan(T(:, col)), col)); % Time to the best solution
    all_rows = true(size(T, 1), 1);
    stats.time_best_c = best(18, solved_c);
    stats.time_best_d = best(19, all_rows);
    stats.time_best_ga = best(20, all_rows);
    stats.time_par_best_cbba = best(21, all_rows);
    stats.time_best_cbba = best(22, all_rows);
    stats.messages_cbba = stat4(T(:, 17));
    stats.err_d = stat4(ok(err_d));
    stats.err_alg1 = stat4(ok(err_alg1));
    stats.err_ga = stat4(ok(err_ga));
    stats.err_cbba = stat4(ok(err_cbba));
    stats.converged_cbba = nnz(T(:, 15) == 1);
    stats.feas_cbba = nnz(T(:, 16) == 1);
    [stats.vs_d.better, stats.vs_d.equal, stats.vs_d.worse] = compare(T(:, 11), T(:, 5));
    [stats.vs_ga.better, stats.vs_ga.equal, stats.vs_ga.worse] = compare(T(:, 11), T(:, 9));
    stats.names = names;
    stats.table = T;

    fprintf("\n> Scenarios: %d (central MILP solved in %d; the cost errors are computed for these)\n", stats.num_tests, stats.num_central);
    fprintf("                                        %10s %10s %10s %10s\n", "Avrg.", "Med.", "Max.", "Min.");
    fprintf("> Comp. time central MILP (9) [s]       %10.3f %10.3f %10.3f %10.3f\n", stats.time_c);
    fprintf("> Comp. time PReDuaL [s]                %10.3f %10.3f %10.3f %10.3f   (Section VII-A)\n", stats.time_d);
    fprintf("> Comp. time Alg. 1 [s]                 %10.3f %10.3f %10.3f %10.3f   (without recovery)\n", stats.time_alg1);
    if has_ga
        fprintf("> Comp. time GA [s]                     %10.3f %10.3f %10.3f %10.3f\n", stats.time_ga);
    end
    fprintf("> Comp. time CBBA, distributed [s]      %10.4f %10.4f %10.4f %10.4f\n", stats.time_par_cbba);
    fprintf("> Comp. time CBBA, serial [s]           %10.4f %10.4f %10.4f %10.4f\n", stats.time_cbba);
    fprintf("> Rounds of CBBA                        %10.1f %10.1f %10d %10d\n", stats.rounds_cbba);
    if ~all(isnan([stats.time_best_c, stats.time_best_d, stats.time_best_ga, stats.time_par_best_cbba]))
        fprintf("> Time to the best solution [s] (NaN: not recorded in the results):\n");
        fprintf(">   Central MILP (9)                    %10.3f %10.3f %10.3f %10.3f\n", stats.time_best_c);
        fprintf(">   PReDuaL (as its comp. time)         %10.3f %10.3f %10.3f %10.3f\n", stats.time_best_d);
        if has_ga
            fprintf(">   GA                                  %10.3f %10.3f %10.3f %10.3f\n", stats.time_best_ga);
        end
        fprintf(">   CBBA, distributed                   %10.4f %10.4f %10.4f %10.4f\n", stats.time_par_best_cbba);
        fprintf(">   CBBA, serial                        %10.4f %10.4f %10.4f %10.4f\n", stats.time_best_cbba);
    end
    fprintf("> Cost error PReDuaL (%%)                %10.2f %10.2f %10.2f %10.2f\n", stats.err_d);
    fprintf("> Cost error Alg. 1 (%%)                 %10.2f %10.2f %10.2f %10.2f   (without recovery)\n", stats.err_alg1);
    if has_ga
        fprintf("> Cost error GA (%%)                     %10.2f %10.2f %10.2f %10.2f\n", stats.err_ga);
    end
    fprintf("> Cost error CBBA (%%)                   %10.2f %10.2f %10.2f %10.2f\n", stats.err_cbba);
    fprintf("> CBBA: converged in %d and feasible in %d of %d scenarios\n", stats.converged_cbba, stats.feas_cbba, stats.num_tests);
    fprintf("> CBBA compared with PReDuaL: better in %d, equal in %d, worse in %d\n", stats.vs_d.better, stats.vs_d.equal, stats.vs_d.worse);
    if has_ga
        fprintf("> CBBA compared with the GA: better in %d, equal in %d, worse in %d\n", stats.vs_ga.better, stats.vs_ga.equal, stats.vs_ga.worse);
    end

end

function [better, equal, worse] = compare(V, V_ref)
    % Number of scenarios where V is better than, equal to and worse than V_ref (relative tolerance 1e-6)
    has = isfinite(V) & isfinite(V_ref);
    tol = 1e-6*max(1, abs(V_ref));
    better = nnz(has & V < V_ref - tol);
    equal = nnz(has & abs(V - V_ref) <= tol);
    worse = nnz(has & V > V_ref + tol);
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
