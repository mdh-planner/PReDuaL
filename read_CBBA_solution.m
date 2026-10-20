%% read_CBBA_solution() - Reads the solution file written by the CBBA executable (see cbba_MRTA)
%
% INPUTS:
% @file: Name of the solution file
%
% OUTPUTS:
% @sol: Structure with fields
%   .name: Name of the instance
%   .feasible: true if every task is assigned to exactly one robot that can perform it (routes are simple paths)
%   .converged: true if CBBA converged (all the robots agree on the winners and no bundle changes)
%   .cost: Cost (8) of the solution, computed by CBBA
%   .time: Serial computation time [s] (all the robots one after the other, wall-clock)
%   .time_par: Distributed computation time [s]: sum over the rounds of the largest time of a robot
%   .time_best, .time_par_best: Serial and distributed time [s] at which the final solution was known (end of
%                               the last round that changed the bundles or winners; NaN if not in the file)
%   .rounds: Number of communication rounds
%   .rounds_best: Last round that changed the bundles or winners
%   .messages: Number of messages sent by the robots
%   .network: Communication network ("full", "ring", "line", "star" or "random")
%   .diameter: Diameter of the communication network
%   .routes: Cell array with the route of each robot (row vector of location indices, starting at its location)
%   .tasks: Cell array with the tasks of each robot (indices of the tasks of the instance file, in the order
%           in which the robot performs them)
%

function sol = read_CBBA_solution(file)
    tok = regexp(fileread(file), '\S+', 'match');
    pos = 1;

    pos = expect(tok, pos, 'MRTA_CBBA_SOLUTION', file);
    [version, pos] = next_number(tok, pos, 'version', file);
    if version ~= 1
        error("read_CBBA_solution:format", "Unsupported version of the solution file %s", file);
    end
    pos = expect(tok, pos, 'name', file);
    [sol.name, pos] = next_token(tok, pos, 'name', file);
    % Numeric fields "key value" until the network (optional fields are NaN if absent)
    fields = {'feasible', 'converged', 'cost', 'time', 'time_par', 'rounds', 'messages'};
    optional = {'time_best', 'time_par_best', 'rounds_best'};
    for q = 1:numel(optional); sol.(optional{q}) = NaN; end
    while pos <= numel(tok) && ~strcmp(tok{pos}, 'network')
        key = tok{pos};
        if ~any(strcmp(key, [fields, optional]))
            error("read_CBBA_solution:format", "Unknown field '%s' in the file %s", key, file);
        end
        [sol.(key), pos] = next_number(tok, pos + 1, key, file);
    end
    missing = fields(~isfield(sol, fields));
    if ~isempty(missing)
        error("read_CBBA_solution:format", "Missing field '%s' in the file %s", missing{1}, file);
    end
    sol.feasible = sol.feasible == 1;
    sol.converged = sol.converged == 1;
    pos = expect(tok, pos, 'network', file);
    [sol.network, pos] = next_token(tok, pos, 'network', file);
    [sol.diameter, pos] = next_number(tok, pos, 'diameter', file);
    pos = expect(tok, pos, 'robots', file);
    [m, pos] = next_number(tok, pos, 'number of robots', file);
    sol.routes = cell(1, m);
    sol.tasks = cell(1, m);
    for i = 1:m
        [sol.routes{i}, pos] = next_list(tok, pos, 'route', file);
        [sol.tasks{i}, pos] = next_list(tok, pos, 'tasks', file);
    end
    expect(tok, pos, 'end', file);

end

function [t, pos] = next_token(tok, pos, what, file)
    % Returns the next token of the file
    if pos > numel(tok)
        error("read_CBBA_solution:format", "Unexpected end of the file %s (expected %s)", file, what);
    end
    t = tok{pos};
    pos = pos + 1;
end

function pos = expect(tok, pos, keyword, file)
    % Reads the next token and checks that it is @keyword
    [t, pos] = next_token(tok, pos, keyword, file);
    if ~strcmp(t, keyword)
        error("read_CBBA_solution:format", "Expected '%s' but found '%s' in the file %s", keyword, t, file);
    end
end

function [v, pos] = next_number(tok, pos, what, file)
    % Reads the next token as a number
    [t, pos] = next_token(tok, pos, what, file);
    v = str2double(t);
    if isnan(v)
        error("read_CBBA_solution:format", "Invalid %s '%s' in the file %s", what, t, file);
    end
end

function [list, pos] = next_list(tok, pos, what, file)
    % Reads a length L followed by L numbers (row vector)
    [L, pos] = next_number(tok, pos, [what, ' length'], file);
    list = zeros(1, L);
    for q = 1:L
        [list(q), pos] = next_number(tok, pos, what, file);
    end
end
