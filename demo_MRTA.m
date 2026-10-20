%% demo_MRTA - Solves a small random MRTA problem with PReDuaL and as the central MILP (9)
%
% Uses MOSEK if it is installed, and MATLAB's intlinprog (Optimization Toolbox) otherwise.
% To force a solver use, e.g., heur_distLag_MRTA(..., solver=struct('solver', "intlinprog")).
clear; clc;

%% Random 4 x 3 grid with m = 4 agents and k = 3 task-types
n_rows = 4;
n_cols = 3;
opt_gen = struct('Evt_sparse', 0.65, 'Eat_sparse', 0.7, 'max_tv', 4, 'min_tv', 1); % Generation options
[Gv, Evt, Eat, names, param] = gen_grid_MRTA(n_rows, n_cols, n_rows, n_cols, opt_gen, 12);
[t_loc, a_loc] = gen_rand_state(Evt, names, opt_gen, 155);

prob = MRTA(Gv, Evt, Eat, names, 'param', param);
prob.normalize_costs(); % Cost normalization (Section VII)
fprintf("Using MILP solver: %s\n", MIQPholder.default_solver());

%% PReDuaL: Algorithm 1 followed by the primal recovery
recovery = true;        % false: Algorithm 1 alone
route_search = true;    % Route search in the primal recovery
use_parallel = false;   % true: the MILP problems of the robots are solved in parallel (Parallel Computing Toolbox)
[z_d, V_d, info_d] = heur_distLag_MRTA(prob, a_loc, t_loc, verbose=true, recovery=recovery, ...
    route_search=route_search, parallel=use_parallel);
fprintf("\nAlgorithm 1: exit flag %d, %d iterations, cost %.4f, time %.3f s\n", ...
    info_d.e_flag_alg1, info_d.k, info_d.V_alg1, info_d.run_time_par_total_alg1);
fprintf("PReDuaL: exit flag %d, cost %.4f, time %.3f s\n", ...
    info_d.e_flag, V_d, info_d.run_time_par_total);
if isfinite(info_d.LB)
    fprintf("Lower bound of the optimal cost: %.4f (optimality gap %.2f%%)\n", info_d.LB, 100*info_d.gap);
end
fprintf("Measured wall-clock time: %.3f s (parallel: %d)\n\n", info_d.run_time_wall_total, info_d.parallel);
prob.state2tex();
prob.res2tex(z_d);

%% Central MILP problem (9)
prob.build_MILP(a_loc, t_loc);
[z_c, V_c, info_c] = prob.miqp.solve();
fprintf("Central MILP: exit flag %d, cost %.4f, time %.3f s\n", info_c.e_flag, V_c, info_c.run_time);
fprintf("Cost error of Algorithm 1: %.2f%%\n", 100*(info_d.V_alg1 - V_c)/V_c);
fprintf("Cost error of PReDuaL: %.2f%%\n", 100*(V_d - V_c)/V_c);

%% Plot the path and tasks of each agent (one figure per agent)
prob.plot_result(z_d, 'plot_weights', false);
