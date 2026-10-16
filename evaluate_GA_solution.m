%% evaluate_GA_solution() - Checks a solution of the GA against the central MILP (9) and computes its cost (8)
%
% Builds the solution z of the central MILP (9) given by the routes and tasks of the GA (Agent.route_to_z),
% checks that it satisfies all the constraints of (9) (MIQPholder.is_feasible), and computes its cost c'*z.
% Note that it builds the central MILP of @prob (prob.build_MILP).
%
% INPUTS:
% @prob: Instance of MRTA (exported with export_GA_instance)
% @a_loc, @t_loc: Locations of the agents and pending tasks
% @task_v, @task_t: Location and type of each task of the instance file (outputs of export_GA_instance)
% @sol: Solution of the GA (output of read_GA_solution)
%
% OUTPUTS:
% @ev: Structure with fields
%   .valid: true if the routes and tasks of the GA are consistent: each route is a simple path of the
%           location graph that starts at the location of the robot, each task is performed once, by a
%           robot that can perform it and visits its location
%   .feasible: true if z satisfies the constraints of (9)
%   .viol: Largest constraint violation of z (see MIQPholder.is_feasible)
%   .V: Cost (8) of z (Inf if it is not feasible)
%   .z: Solution of (9) ([] if not valid)
%   .msg: Description of the first problem found ('' if none)
%

function ev = evaluate_GA_solution(prob, a_loc, t_loc, task_v, task_t, sol)
    ev = struct('valid', false, 'feasible', false, 'viol', inf, 'V', inf, 'z', [], 'msg', '');
    prob.build_MILP(a_loc, t_loc, 'sparse', true, 'condense', true);
    m = prob.m;
    n_t = length(task_v);

    if numel(sol.routes) ~= m || numel(sol.tasks) ~= m
        ev.msg = sprintf('the solution has %d robots instead of %d', numel(sol.routes), m);
        return
    end
    all_tasks = [sol.tasks{:}];
    if numel(all_tasks) ~= n_t || ~isequal(sort(all_tasks(:))', 1:n_t)
        ev.msg = 'not every task is performed exactly once';
        return
    end

    z = cell(m, 1);
    for i = 1:m
        ag = prob.agent{i};
        route = sol.routes{i};
        tasks = sol.tasks{i};
        if isempty(route) || any(route < 1 | route > prob.n | route ~= round(route)) || route(1) ~= ag.loc
            ev.msg = sprintf('invalid route of robot %d', i);
            return
        end
        if numel(unique(route)) ~= numel(route)
            ev.msg = sprintf('the route of robot %d visits a location more than once', i);
            return
        end
        if length(route) > 1 && any(findedge(prob.Gv, route(1:end-1), route(2:end)) == 0)
            ev.msg = sprintf('the route of robot %d uses an edge that is not in the location graph', i);
            return
        end
        tasks_on = zeros(1, numel(tasks));
        for q = 1:numel(tasks)
            v = task_v(tasks(q));
            t = task_t(tasks(q));
            tasks_on(q) = ag.idx_Evt(v, t);
            if tasks_on(q) == 0
                ev.msg = sprintf('robot %d cannot perform task %d', i, tasks(q));
                return
            end
            if ~any(route == v)
                ev.msg = sprintf('robot %d does not visit the location of task %d', i, tasks(q));
                return
            end
        end
        z{i} = ag.route_to_z(route, tasks_on);
    end

    ev.valid = true;
    ev.z = cat(1, z{:});
    [ev.feasible, ev.viol] = prob.miqp.is_feasible(ev.z);
    if ev.feasible
        ev.V = full(prob.miqp.primal.c'*ev.z);
    else
        ev.msg = sprintf('the solution does not satisfy the constraints of (9) (violation %g)', ev.viol);
    end

end
