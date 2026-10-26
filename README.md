# PReDuaL: a primal-recovering residual-driven dual Lagrangian method for MRTA

This MATLAB package implements PReDuaL, a method for multi-robot task allocation (MRTA) problems in which pending tasks are located at the nodes of a directed location graph and each robot must be assigned tasks and a route to perform them.

The problem is posed as a mixed-integer linear program (MILP). A dual Lagrangian method decomposes it into one small MILP per robot, which can be solved in parallel, and a heuristic update of the dual variables driven by the residual of the coupling constraints steers the iterates towards feasible task allocations (Algorithms 1 and 2 of the article). A primal recovery step then combines the routes computed by the robots into a feasible solution: column generation over the pool of routes, which also gives a lower bound of the optimal cost, followed by route selection, route search and re-routing (Algorithm 3).

## Requirements

- MATLAB R2021a or newer.
- A MILP solver: [MOSEK](https://docs.mosek.com/latest/toolbox/install-interface.html) (used in the article, tested with MOSEK 10.2), or MATLAB's `intlinprog` (Optimization Toolbox), which is used if MOSEK is not found.
- Optional: the Parallel Computing Toolbox, to solve the MILP problems of the robots in parallel (`parallel=true`).
- The comparison with the genetic algorithm (GA) and the consensus-based bundle algorithm (CBBA) uses the executables `bin/ga_mrta.exe` and `bin/cbba_mrta.exe`, which run on 64-bit Windows. On other systems, set `use_ga = false` and `use_cbba = false` in `main_run`.

No installation is needed beyond adding the folder to the MATLAB path.

## Quick start

`demo_MRTA` generates a small random problem, solves it with PReDuaL and as a central MILP, and plots the route and tasks of each robot:

```matlab
opt = struct('Evt_sparse', 0.65, 'Eat_sparse', 0.7, 'max_tv', 4, 'min_tv', 1);
[Gv, Evt, Eat, names, param] = gen_grid_MRTA(4, 3, 4, 3, opt, 12); % 4 x 3 grid, 4 robots, 3 task-types
[t_loc, a_loc] = gen_rand_state(Evt, names, opt, 155);           % Robot locations and pending tasks
prob = MRTA(Gv, Evt, Eat, names, 'param', param);
[z, V, info] = heur_distLag_MRTA(prob, a_loc, t_loc, verbose=true);
prob.res2tex(z);
```

## Reproducing the results of the article

Run `main_run`. For each grid from 8 x 7 to 13 x 12, it solves 10 random scenarios with PReDuaL, with MOSEK (time limit of 3600 s), with the GA given the same computation time as PReDuaL, and with CBBA. The results are saved in `results/` and summarized by

```matlab
analyze_rand_batch(n);                     % PReDuaL and the central MILP, n = number of locations, e.g., 56 for the 8 x 7 grid
analyze_GA_batch(n, time_budget="paper");  % GA vs PReDuaL and the central MILP
analyze_CBBA_batch(n);                     % CBBA vs PReDuaL, the GA and the central MILP
```

An interrupted run can be resumed by running `main_run` again (`skip_done=true` skips the scenarios already solved).

## Correspondence with the article

| Article | Code |
| --- | --- |
| Constraints (1), (4)-(7) and cost (8) of each robot (set Z_i) | `Agent.build_MILP` |
| Constraints (7) | `Agent.gen_valid_ineq_constraints` (option `tight`) |
| Coupling constraints (3) / (10b) | `MRTA.gen_single_task_eq_constraints` |
| Central MILP (9) | `MRTA.build_MILP` |
| Algorithm 1 | `heur_distLag_MRTA` |
| Algorithm 2 | `heuristic_residual` (in `heur_distLag_MRTA.m`) |
| Times Gamma_i(t_j) | `Agent.task_times` |
| Primal recovery (Section VI, Algorithm 3) | `primal_recovery_MRTA`, `RoutePool`, `RouteGraph`, `RouteEvaluator` |
| Route selection (18) and its LP relaxation (19) | `RoutePool.selection_problem` |
| Column generation, pricing problems (20) and lower bound (21) | `column_generation`, `price_all` (in `primal_recovery_MRTA.m`), `RouteEvaluator.price_agents` |
| Route search, routes (a)-(d) (Section VI-C) | `RouteGraph` |
| Re-routing with Held-Karp (Section VI-C) | `RouteEvaluator.route_dp` |
| Implementation details (Section VII) | `MRTA.normalize_costs`, `Agent.prize_bounds`, `Agent.warm_start_z` |
| Random grid problems (Section VII-A) | `gen_rand_grid_instance`, `gen_grid_MRTA`, `gen_rand_state`, `gen_param_grid` |
| Computation time (Sections VI-D and VII-A) | `info.run_time_par_total` returned by `heur_distLag_MRTA` |
| Genetic algorithm (Section VII-C) | `bin/ga_mrta.exe`, called by `ga_MRTA` (`run_GA_batch`, `analyze_GA_batch`) |
| CBBA (Section VII-D) | `bin/cbba_mrta.exe`, called by `cbba_MRTA` (`run_CBBA_batch`, `analyze_CBBA_batch`) |

## Main options of `heur_distLag_MRTA`

The defaults are the settings of the article (Section VII-A). They can be changed with name/value arguments, e.g., `heur_distLag_MRTA(prob, a_loc, t_loc, k_f=10, inner_max_time=20)`:

- `recovery`: apply the primal recovery after Algorithm 1 (default `true`).
- `column_generation`: start the primal recovery with column generation on the routes of Algorithm 1 (default `true`). It gives a lower bound of the optimal cost (`info.LB`, `info.gap`).
- `route_search`: apply the route search in the primal recovery (default `true`).
- `recovery_rounds`: number of rounds N_rec of the primal recovery (default 1).
- `inner_max_time`: time limit of the MILP problems of the robots in Algorithm 1 (default 10 s).
- `cg_max_time`: time limit of the pricing problems in the column generation (default 2 s). When the column generation would stop, it switches to `inner_max_time` and continues until it stops again.
- `cg_stall`, `cg_max_iter`, `cg_smoothing`, `cg_tol`: k_cg, N_cg, zeta and epsilon_cg of the column generation (defaults 3, 20, 0.5 and 1e-4).
- `tight`, `prize_bounds`, `warm_start`, `pool_stall`: constraints (7) and the implementation details of Sections V and VII (default `true`, `true`, `true` and k_s = 3).
- `parallel`: solve the MILP problems of the robots in parallel with `parfor` (default `false`).
- `solver`: options of the MILP solver, e.g., `solver=struct('solver', "intlinprog")`.

See the header of `heur_distLag_MRTA.m` for the complete list of options and outputs.

## Computation time

As in the article, the robots are assumed to solve their MILP problems in parallel, so the computation time is `overhead + sum_k max_i (time of the MILP problem of robot i at iteration k)`, and the primal recovery is timed in the same way (`info.run_time_par_total`). By default, the MILP problems are solved one after the other, so the measured wall-clock time (`info.run_time_wall_total`) is larger. With `parallel=true`, they are solved in parallel and the wall-clock time approaches the computation time of the article on a computer with at least as many cores as robots.

`run_rand_grid` also saves the time at which each method found the solution it returns: the central MILP (`res.time_best_c`, from the log of MOSEK or an output function of intlinprog; `res.hist_c` has the time and cost of each new best solution) and PReDuaL (`res.time_best_d`, computed as `res.time_d_par_total` up to that point; `res.hist_d` has the new best solutions of Algorithm 1 and the primal recovery). The GA (`res_ga.time_best`) and CBBA (`res_cbba.time_par_best`, `res_cbba.time_best`) record it too, and the analysis functions report it.

> [!NOTE]
> This repository is the replication package of the article "PReDuaL: A Primal-Recovering Residual-Driven Dual Lagrangian Method for Multi-Robot Task Allocation" (under revision).
