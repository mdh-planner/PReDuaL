%% export_GA_instance() - Writes an MRTA problem to the input file of the GA and CBBA executables (see ga_MRTA)
%
% The file contains the directed location graph (edges and their weights w_{i,j}, with the costs of @prob,
% i.e., normalized if prob.normalize_costs() was called), the location and capabilities (task types it can
% perform) of each robot, the pending tasks (location, type and cost rho), and the costs w_{i,i} (cost_stop)
% and epsilon_tau (cost_times) of the cost function (8).
%
% The pending tasks are the task-location pairs of the condensed problem (see MRTA.condense_Evt), in the
% order of the coupling constraints (3) of MRTA.gen_single_task_eq_constraints (by task type, then location).
%
% INPUTS:
% @prob: Instance of MRTA
% @a_loc: Cell array with the location (name) of each agent
% @t_loc: Cell array with the locations (names) of the pending tasks of each task type
% @file: Name of the file
% @name: (optional) Name of the instance (default "instance")
%
% OUTPUTS:
% @task_v, @task_t: Location and type of each pending task (task j of the file)
%

function [task_v, task_t] = export_GA_instance(prob, a_loc, t_loc, file, name)
    if nargin < 5; name = "instance"; end
    name = regexprep(char(name), '\s', '_');
    if isempty(name); name = 'instance'; end
    if length(a_loc) ~= prob.m
        error("export_GA_instance:a_loc", "The location of all the agents must be specified");
    end

    % Pending tasks (as MRTA.condense_Evt, without modifying @prob)
    Evt_p = zeros(prob.n, prob.k);
    for t = 1:prob.k
        t_loc_t = t_loc{t};
        for j = 1:length(t_loc_t)
            v = prob.v_pos(t_loc_t{j});
            if isempty(v)
                error("export_GA_instance:t_loc", "Unknown location %s of a task of type %d", char(t_loc_t{j}), t);
            end
            if prob.Evt_full(v, t)
                Evt_p(v, t) = 1;
            end
        end
    end
    [task_v, task_t] = find(Evt_p); % Ordered by task type, then location (as the coupling constraints (3))
    task_v = reshape(task_v, [], 1);
    task_t = reshape(task_t, [], 1);

    % Location graph
    [e_s, e_t] = findedge(prob.Gv);
    e_w = prob.Gv.Edges.Weight;

    fid = fopen(file, 'w');
    if fid < 0
        error("export_GA_instance:file", "Cannot write the file %s", file);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'MRTA_GA_INSTANCE 1\n');
    fprintf(fid, 'name %s\n', name);
    fprintf(fid, 'locations %d\n', prob.n);
    fprintf(fid, 'types %d\n', prob.k);
    fprintf(fid, 'edges %d\n', length(e_w));
    fprintf(fid, '%d %d %.17g\n', [e_s(:)'; e_t(:)'; e_w(:)']);
    fprintf(fid, 'robots %d\n', prob.m);
    for i = 1:prob.m
        loc = prob.v_pos(a_loc{i});
        if isempty(loc)
            error("export_GA_instance:a_loc", "Unknown location of agent %d", i);
        end
        fprintf(fid, '%d', loc);
        fprintf(fid, ' %d', double(prob.Eat(i, :) ~= 0));
        fprintf(fid, '\n');
    end
    fprintf(fid, 'tasks %d\n', length(task_v));
    for j = 1:length(task_v)
        fprintf(fid, '%d %d %.17g\n', task_v(j), task_t(j), prob.param.cost_task(task_t(j)));
    end
    fprintf(fid, 'cost_stop %.17g\n', prob.param.cost_stop);
    fprintf(fid, 'cost_times %.17g\n', prob.param.cost_times);
    fprintf(fid, 'end\n');

end
