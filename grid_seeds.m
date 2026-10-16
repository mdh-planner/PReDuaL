%% grid_seeds() - State seeds of the random scenarios of a grid
%
% seed_state = grid_seeds(n, num_tests)
%
% Returns num_tests different state seeds in [1, 1000] for the scenarios of a grid with n locations
% (used by run_grid_batch and run_GA_batch). A repeated seed, which would give the same scenario twice,
% is replaced by a new random draw.
%
% INPUTS:
% @n: Number of locations of the grid (n_rows * n_cols), used as RNG seed
% @num_tests: Number of scenarios
%
% OUTPUTS:
% @seed_state: 1 x num_tests vector of different seeds
%

function seed_state = grid_seeds(n, num_tests)
    rng(n);
    seed_state = randi([1, 1000], [1, num_tests]);
    for i = 2:num_tests
        while any(seed_state(1:i-1) == seed_state(i))
            seed_state(i) = randi([1, 1000]);
        end
    end
end
