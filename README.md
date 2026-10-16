# PReDuaL: a primal-recovering residual-driven dual Lagrangian method for MRTA

This MATLAB package implements PReDuaL, a method for multi-robot task allocation (MRTA) problems in which pending tasks are located at the nodes of a directed location graph and each robot must be assigned tasks and a route to perform them.

The problem is posed as a mixed-integer linear program (MILP). A dual Lagrangian method decomposes it into one small MILP per robot, which can be solved in parallel, and a heuristic update of the dual variables driven by the residual of the coupling constraints steers the iterates towards feasible task allocations (Algorithms 1 and 2 of the article). A primal recovery step then combines the routes computed by the robots during the iterations into a feasible solution with a cost close to the optimal one (Algorithm 3).

## Requirements

- MATLAB R2021a or newer.
- A MILP solver: [MOSEK](https://docs.mosek.com/latest/toolbox/install-interface.html) (used in the article, tested with MOSEK 10.2), or MATLAB's `intlinprog` (Optimization Toolbox), which is used if MOSEK is not found.
- Optional: the Parallel Computing Toolbox, to solve the MILP problems of the robots in parallel (`parallel=true`).
- The comparison with the genetic algorithm (GA) uses the executable `bin/ga_mrta.exe`, which runs on 64-bit Windows. On other systems, set `use_ga = false` in `main_run`.

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

Run `main_run`. For each grid from 8 x 7 to 12 x 11, it solves 10 random scenarios with PReDuaL and with MOSEK (time limit of 3600 s), and with the GA given the same computation time as PReDuaL. The results are saved in `results/` and summarized as in Tables I and II by

```matlab
analyze_rand_batch(n);                     % Table I, n = number of locations, e.g., 56 for the 8 x 7 grid
analyze_GA_batch(n, time_budget="paper");  % Table II
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
| Preconditioning and solution of the inner problems (Section VII-A) | `MRTA.normalize_costs`, `Agent.prize_bounds`, `Agent.warm_start_z` |
| Random grid problems (Section VII-B) | `gen_rand_grid_instance`, `gen_grid_MRTA`, `gen_rand_state`, `gen_param_grid` |
| Computation time (Section VII-B) | `info.run_time_par_total` returned by `heur_distLag_MRTA` |
| Genetic algorithm (Section VII-D) | `bin/ga_mrta.exe`, called by `ga_MRTA` (`run_GA_batch`, `analyze_GA_batch`) |

## Main options of `heur_distLag_MRTA`

The default hyperparameters of Algorithms 1 and 2 are those of Section VII-B. They, and the following switches, can be changed with name/value arguments, e.g., `heur_distLag_MRTA(prob, a_loc, t_loc, k_f=10, inner_max_time=10)`:

- `recovery`: apply the primal recovery after Algorithm 1 (default `true`).
- `light_recovery`: only route selection, route search and re-routing (default `true`). With `false`, the local search is also applied, as in the article.
- `recovery_rounds`: number of rounds N_rec of the primal recovery (default 3; 1 in the article).
- `column_generation`: with `light_recovery=false`, apply column generation after the primal recovery (default `true`; `false` in the article). It finds new routes and gives a lower bound of the optimal cost (`info.LB`, `info.gap`).
- `pair_reopt`: re-optimize the tasks of each pair of robots in the primal recovery (default `false`).
- `inner_max_time`: time limit of the MILP problems of the robots in Algorithm 1 (default `inf`; 10 s in the article).
- `tight`, `prize_bounds`, `warm_start`, `pool_stall`: formulation (7) and the implementation details of Sections V and VII-A (default `true`, `true`, `true` and 3). Use `tight=false, prize_bounds=false, warm_start=false, pool_stall=0` for the formulation and method of the original submission.
- `parallel`: solve the MILP problems of the robots in parallel with `parfor` (default `false`).
- `solver`: options of the MILP solver, e.g., `solver=struct('solver', "intlinprog")`.

`main_run` uses the settings of the article. See the header of `heur_distLag_MRTA.m` for the complete list of options and outputs.

## Computation time

As in the article, the robots are assumed to solve their MILP problems in parallel, so the computation time is `overhead + sum_k max_i (time of the MILP problem of robot i at iteration k)`, and the primal recovery is timed in the same way (`info.run_time_par_total`). By default, the MILP problems are solved one after the other, so the measured wall-clock time (`info.run_time_wall_total`) is larger. With `parallel=true`, they are solved in parallel and the wall-clock time approaches the computation time of the article on a computer with at least as many cores as robots.

> [!NOTE]
> This repository is the replication package of the article "PReDuaL: A Primal-Recovering Residual-Driven Dual Lagrangian Method for Multi-Robot Task Allocation" by B. Miloradović, P. Krupa, M. Frasheri and C. Trubiani (under revision).
