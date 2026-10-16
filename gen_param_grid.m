%% gen_param_grid - Returns the parameters of the random problems of Section VII-B
%
% INPUTS:
%   @n_rows: int, number of rows of the grid
%   @n_cols: int, number of columns of the grid
%   @m: int, number of agents
%   @k: int, number of task-types
%
% OUTPUTS:
%   @param: struct, containing the parameters for the random problem generation
%

function param = gen_param_grid(n_rows, n_cols, m, k)

    param.row_lim = [n_rows, n_rows]; % Force dimension of grid
    param.col_lim = [n_cols, n_cols]; % Force dimension of grid
    param.m_lim = [m, m]; % Force number of agents
    param.k_lim = [k, k]; % Force number of task-types
    param.Evt_sparse_lim = [0.5, 0.8]; % Limits of P_V (task-types of the locations)
    param.Eat_sparse_lim = [0.5, 0.9]; % Limits of P_S (task-types of the robots)
    param.max_tv_lim = [3, 6]; % Limits of n_max (maximum number of tasks of each type)
    param.min_tv_lim = [0, 2]; % Limits of n_min (minimum number of tasks of each type)

end
