%% ga_MRTA() - Solves an MRTA problem with the genetic algorithm of Section VII-D
%
% [z, V, info] = ga_MRTA(prob, a_loc, t_loc, Name=Value)
%
% Writes the problem to a text file (export_GA_instance), solves it with the GA executable bin/ga_mrta.exe,
% reads its solution (read_GA_solution) and checks it against the central MILP (9)
% (evaluate_GA_solution). The GA solves the same problem as PReDuaL: same cost (8), each pending task
% is performed once by a robot that can perform it, and the route of each robot is a simple path.
%
% INPUTS:
% @prob: Instance of MRTA (with the same costs as for the other solvers, e.g., after prob.normalize_costs())
% @a_loc: Cell array with the location (name) of each agent
% @t_loc: Cell array with the locations (names) of the pending tasks of each task type
%
% NAME/VALUE INPUTS:
% 'exe': Path of the GA executable (default: bin/ga_mrta.exe of this package)
% 'seed': Seed of the random generator of the GA (default 1; 0 for a random seed)
% 'pop': Population size (default 200)
% 'stop': The GA stops after 'stop' generations without improvement (default 200)
% 'maxgen': Maximum number of generations (default 100000)
% 'time': Time limit of the GA [s] (default inf)
% 'restarts': Time budget mode (requires 'time'): when a run of the GA stops ('stop' generations without
%             improvement), a new run is started (new initial population) while there is time left, and the
%             best solution of all the runs is returned (default false: a single run)
% 'args': Additional options of the GA (text): "--mutation P", "--crossover P" and "--elite P" (defaults
%         0.15, 0.7 and 0.05), e.g., "--mutation 0.2"
% 'folder': Folder for the instance and solution files (default: tempdir)
% 'name': Name of the instance (default "instance")
% 'keep_files': Keep the instance and solution files (default false)
% 'verbose': Print the output of the GA (default false)
%
% OUTPUTS:
% @z: Solution of the central MILP (9) given by the GA ([] if it is not valid)
% @V: Cost (8) of z (Inf if it is not feasible)
% @info: Structure with fields
%   .e_flag: 1 if the GA found a feasible solution of (9), -1 otherwise
%   .feasible: true if z is a feasible solution of (9)
%   .time: Computation time of the GA [s] (wall-clock time measured by the executable: initial population
%          and generations, without reading/writing files)
%   .wall_time: Wall-clock time of the call to the executable [s] (including starting it and the files)
%   .generations: Number of generations of the GA (sum over its runs)
%   .runs: Number of runs of the GA (more than 1 only with 'restarts')
%   .cost_ga: Cost (8) computed by the GA (equal to V up to rounding errors)
%   .routes, .tasks: Route (locations) and tasks of each robot (see read_GA_solution)
%   .task_v, .task_t: Location and type of each task (see export_GA_instance)
%   .viol, .msg: See evaluate_GA_solution
%   .files: Instance and solution files (deleted unless 'keep_files' is true)
%

function [z, V, info] = ga_MRTA(prob, a_loc, t_loc, opt)
    arguments
        prob MRTA
        a_loc
        t_loc
        opt.exe {mustBeText} = ""
        opt.seed (1, 1) {mustBeNonnegative, mustBeInteger} = 1
        opt.pop (1, 1) {mustBePositive, mustBeInteger} = 200
        opt.stop (1, 1) {mustBePositive, mustBeInteger} = 200
        opt.maxgen (1, 1) {mustBePositive, mustBeInteger} = 100000
        opt.time (1, 1) double {mustBePositive} = inf
        opt.restarts (1, 1) {mustBeNumericOrLogical} = false
        opt.args {mustBeText} = ""
        opt.folder {mustBeText} = ""
        opt.name {mustBeText} = "instance"
        opt.keep_files (1, 1) {mustBeNumericOrLogical} = false
        opt.verbose (1, 1) {mustBeNumericOrLogical} = false
    end

    if opt.restarts && ~isfinite(opt.time)
        error("ga_MRTA:restarts", "The option 'restarts' requires a time limit ('time')");
    end
    exe = ga_executable(opt.exe);
    folder = char(opt.folder);
    if isempty(folder); folder = tempdir; end
    if ~isfolder(folder); mkdir(folder); end
    base = tempname(folder); % Unique name
    inst_file = [base, '_instance.txt'];
    sol_file = [base, '_solution.txt'];
    info.files = {inst_file, sol_file};

    % Write the instance
    [task_v, task_t] = export_GA_instance(prob, a_loc, t_loc, inst_file, opt.name);

    % Run the GA
    ga_args = {inst_file, sol_file, '--seed', sprintf('%d', opt.seed), '--pop', sprintf('%d', opt.pop), ...
               '--stop', sprintf('%d', opt.stop), '--maxgen', sprintf('%d', opt.maxgen)};
    if isfinite(opt.time)
        ga_args = [ga_args, {'--time', sprintf('%.17g', opt.time)}];
    end
    if opt.restarts
        ga_args = [ga_args, {'--restarts'}];
    end
    if ~opt.verbose
        ga_args = [ga_args, {'--quiet'}];
    end
    extra = strtrim(char(opt.args));
    t_wall = tic;
    [status, out] = system(command_line(exe, ga_args, extra));
    info.wall_time = toc(t_wall);
    if opt.verbose; fprintf("%s", out); end
    if status ~= 0 || ~isfile(sol_file)
        cleanup_files(info.files, opt.keep_files);
        error("ga_MRTA:ga", "The GA (%s) failed with exit code %d:\n%s", exe, status, out);
    end

    % Read and check the solution
    sol = read_GA_solution(sol_file);
    ev = evaluate_GA_solution(prob, a_loc, t_loc, task_v, task_t, sol);
    cleanup_files(info.files, opt.keep_files);

    z = ev.z;
    V = ev.V;
    info.feasible = ev.feasible;
    info.e_flag = 1;
    if ~ev.feasible; info.e_flag = -1; end
    info.time = sol.time;
    info.generations = sol.generations;
    info.runs = sol.runs;
    info.cost_ga = sol.cost;
    info.feasible_ga = sol.feasible;
    info.routes = sol.routes;
    info.tasks = sol.tasks;
    info.task_v = task_v;
    info.task_t = task_t;
    info.viol = ev.viol;
    info.msg = ev.msg;

    if sol.feasible && ~ev.feasible
        warning("ga_MRTA:check", "The solution of the GA is not a feasible solution of (9): %s", ev.msg);
    elseif ev.feasible && abs(V - sol.cost) > 1e-6*max(1, abs(V))
        warning("ga_MRTA:cost", "The cost computed by the GA (%.10g) differs from the cost (8) of its solution (%.10g)", sol.cost, V);
    end

end

function exe = ga_executable(exe)
    % Path of the GA executable
    exe = char(exe);
    if ~isempty(exe)
        if ~isfile(exe)
            error("ga_MRTA:exe", "The GA executable %s was not found", exe);
        end
        return
    end
    bin = fullfile(fileparts(mfilename('fullpath')), 'bin');
    if ispc
        exe = fullfile(bin, 'ga_mrta.exe');
    else
        exe = fullfile(bin, 'ga_mrta'); % Not included in the package (only the Windows executable is)
    end
    if ~isfile(exe)
        error("ga_MRTA:exe", "The GA executable %s was not found (use the option 'exe' to give its path)", exe);
    end
end

function cmd = command_line(exe, args, extra)
    % Command line to run the executable with the given arguments (paths with spaces are quoted)
    parts = [{exe}, args];
    quoted = false;
    for q = 1:numel(parts)
        if any(isspace(parts{q}))
            parts{q} = ['"', parts{q}, '"'];
            quoted = true;
        end
    end
    cmd = strjoin(parts, ' ');
    if ~isempty(extra)
        cmd = [cmd, ' ', extra];
    end
    cmd = [cmd, ' 2>&1'];
    if ispc && quoted
        cmd = ['"', cmd, '"']; % cmd.exe removes the first and last quotes of the command line
    end
end

function cleanup_files(files, keep)
    % Deletes the instance and solution files (unless @keep is true)
    if keep; return; end
    for q = 1:numel(files)
        if isfile(files{q}); delete(files{q}); end
    end
end
