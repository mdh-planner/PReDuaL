%% gen_rand_grid_instance() - Generates a random grid problem (Section VII-A)
%
% Used by run_rand_grid and run_GA_batch, so that all the methods solve the same problems.
%
% INPUTS:
% @param: struct with the parameters for random generation (see gen_param_grid and run_rand_grid)
% @seed: struct with fields .graph and .state (see run_rand_grid)
%
% OUTPUTS:
% @prob: Instance of MRTA with normalized costs (Section VII)
% @a_loc, @t_loc: Locations of the agents and pending tasks
% @info: Structure with fields n_rows, n_cols (size of the grid) and num_tasks (number of pending tasks)
%

function [prob, a_loc, t_loc, info] = gen_rand_grid_instance(param, seed)

    rng(seed.graph); % Set seed for generating random graph

    % Size of grid, agents and task-types
    n_rows = randi(param.row_lim); % Rows of grid
    n_cols = randi(param.col_lim); % Columns of grid
    m = randi(param.m_lim); % Number of agents
    k = randi(param.k_lim); % Number of task-types

    % Probabilities P_V and P_S of Section VII-A
    opt_MRTA.Evt_sparse = (param.Evt_sparse_lim(2) - param.Evt_sparse_lim(1))*rand(1) + param.Evt_sparse_lim(1);
    opt_MRTA.Eat_sparse = (param.Eat_sparse_lim(2) - param.Eat_sparse_lim(1))*rand(1) + param.Eat_sparse_lim(1);

    % Bounds n_max and n_min of the number of tasks of each task-type
    opt_MRTA.max_tv = randi(param.max_tv_lim);
    opt_MRTA.min_tv = randi(param.min_tv_lim);

    % Generate MRTA problem
    [Gv, Evt, Eat, names, param_MRTA] = gen_grid_MRTA(n_rows, n_cols, m, k, opt_MRTA, seed.graph);
    [t_loc, a_loc] = gen_rand_state(Evt, names, opt_MRTA, seed.state);

    prob = MRTA(Gv, Evt, Eat, names, 'param', param_MRTA);
    prob.normalize_costs(); % Cost normalization (Section VII)

    info.n_rows = n_rows;
    info.n_cols = n_cols;
    info.num_tasks = prob.get_num_tasks(t_loc);

end
