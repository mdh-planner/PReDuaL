%% main_run - Reproduces the results of Table I of the article
%   B. Miloradovic, P. Krupa, M. Frasheri, C. Trubiani, "PReDuaL: A Primal-Recovering Residual-Driven
%   Dual Lagrangian Method for Multi-Robot Task Allocation".
%
% Solves the 10 random scenarios of each grid with PReDuaL, the central MILP, the GA given the same
% computation time as PReDuaL, and CBBA. The results are saved in results/.
% The MILP solver can be selected with, e.g., run_grid_batch(..., solver=struct('solver', "intlinprog")).
clear; clc;

grids = 8:13; % Grids of n_v x m_v = n x (n - 1)

%% Switches
recovery = true;            % Primal recovery after Algorithm 1 (false: Algorithm 1 alone)
route_search = true;        % Route search in the primal recovery
col_gen = true;             % Start the primal recovery with column generation (lower bound of the optimal cost)
rec_rounds = 1;             % Number of rounds N_rec of the primal recovery
cg_stall = 3;               % Rounds k_cg without improvement of the LP relaxation that stop the column generation
cg_max_time = 2;            % Time limit [s] of the pricing MILPs (inner_max_time in the last rounds)
inner_max_time = 10;        % Time limit [s] of the MILP problems of the robots in Algorithm 1 (inf: none)
use_parallel = false;       % Solve the MILP problems of the robots in parallel (Parallel Computing Toolbox)
use_ga = true;              % Also solve the scenarios with the GA (bin/ga_mrta.exe, Windows only)
ga_time_budget = "paper";   % "paper": the GA gets the computation time of PReDuaL of Section VII-A on each
                            % scenario; "wall": the measured wall-clock time of PReDuaL; "none": the GA
                            % stops with its own rule (200 generations without improvement)
use_cbba = true;            % Also solve the scenarios with CBBA (bin/cbba_mrta.exe, Windows only)
cbba_opt = struct();        % Options of cbba_MRTA, e.g., struct('network', "ring")
skip_done = true;           % Skip the scenarios already solved (to resume an interrupted run)

heur = struct('recovery', recovery, 'route_search', route_search, 'column_generation', col_gen, ...
              'recovery_rounds', rec_rounds, 'cg_stall', cg_stall, 'cg_max_time', cg_max_time, ...
              'inner_max_time', inner_max_time, 'parallel', use_parallel);

for r = grids
    c = r - 1;
    n = r*c;
    run_grid_batch(r, c, heur=heur, central_tight=true, central_max_time=3600, ...
                   solver=struct('solver', "mosek"), skip_done=skip_done);
    analyze_rand_batch(n);

    if use_ga
        run_GA_batch(r, c, time_budget=ga_time_budget, skip_done=skip_done);
        analyze_GA_batch(n, time_budget=ga_time_budget);
    end

    if use_cbba
        run_CBBA_batch(r, c, cbba=cbba_opt, skip_done=skip_done);
        if use_ga; ga_budget = ga_time_budget; else; ga_budget = ""; end
        analyze_CBBA_batch(n, ga_budget=ga_budget);
    end
end
