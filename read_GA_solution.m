%% read_GA_solution() - Reads the solution file written by the GA executable (see ga_MRTA)
%
% INPUTS:
% @file: Name of the solution file
%
% OUTPUTS:
% @sol: Structure with fields
%   .name: Name of the instance
%   .feasible: true if the GA found a feasible solution
%   .cost: Cost (8) of the solution, computed by the GA
%   .time: Computation time of the GA [s] (wall-clock time measured by the GA)
%   .generations: Number of generations of the GA (sum over its runs)
%   .runs: Number of runs of the GA (more than 1 only with its option --restarts, see ga_MRTA)
%   .seed: Seed of the random generator of the GA
%   .routes: Cell array with the route of each robot (row vector of location indices, starting at its location)
%   .tasks: Cell array with the tasks of each robot (indices of the tasks of the instance file, in the order
%           in which the robot performs them)
%

function sol = read_GA_solution(file)
    txt = fileread(file);
    tok = regexp(txt, '\S+', 'match');
    pos = 1;

    pos = expect(tok, pos, 'MRTA_GA_SOLUTION', file);
    [version, pos] = next_number(tok, pos, 'version', file);
    if version ~= 1
        error("read_GA_solution:format", "Unsupported version of the solution file %s", file);
    end
    pos = expect(tok, pos, 'name', file);
    [sol.name, pos] = next_token(tok, pos, 'name', file);
    pos = expect(tok, pos, 'feasible', file);
    [feasible, pos] = next_number(tok, pos, 'feasible', file);
    sol.feasible = feasible == 1;
    pos = expect(tok, pos, 'cost', file);
    [sol.cost, pos] = next_number(tok, pos, 'cost', file);
    pos = expect(tok, pos, 'time', file);
    [sol.time, pos] = next_number(tok, pos, 'time', file);
    pos = expect(tok, pos, 'generations', file);
    [sol.generations, pos] = next_number(tok, pos, 'generations', file);
    sol.runs = 1;
    if pos <= numel(tok) && strcmp(tok{pos}, 'runs') % Optional field
        [sol.runs, pos] = next_number(tok, pos + 1, 'runs', file);
    end
    pos = expect(tok, pos, 'seed', file);
    [sol.seed, pos] = next_number(tok, pos, 'seed', file);
    pos = expect(tok, pos, 'robots', file);
    [m, pos] = next_number(tok, pos, 'number of robots', file);
    sol.routes = cell(1, m);
    sol.tasks = cell(1, m);
    for i = 1:m
        [L, pos] = next_number(tok, pos, 'route length', file);
        route = zeros(1, L);
        for q = 1:L
            [route(q), pos] = next_number(tok, pos, 'location', file);
        end
        [nt, pos] = next_number(tok, pos, 'number of tasks', file);
        tasks = zeros(1, nt);
        for q = 1:nt
            [tasks(q), pos] = next_number(tok, pos, 'task', file);
        end
        sol.routes{i} = route;
        sol.tasks{i} = tasks;
    end
    expect(tok, pos, 'end', file);

end

function [t, pos] = next_token(tok, pos, what, file)
    % Returns the next token of the file
    if pos > numel(tok)
        error("read_GA_solution:format", "Unexpected end of the file %s (expected %s)", file, what);
    end
    t = tok{pos};
    pos = pos + 1;
end

function pos = expect(tok, pos, keyword, file)
    % Reads the next token and checks that it is @keyword
    [t, pos] = next_token(tok, pos, keyword, file);
    if ~strcmp(t, keyword)
        error("read_GA_solution:format", "Expected '%s' but found '%s' in the file %s", keyword, t, file);
    end
end

function [v, pos] = next_number(tok, pos, what, file)
    % Reads the next token as a number
    [t, pos] = next_token(tok, pos, what, file);
    v = str2double(t);
    if isnan(v)
        error("read_GA_solution:format", "Invalid %s '%s' in the file %s", what, t, file);
    end
end
