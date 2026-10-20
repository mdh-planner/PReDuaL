%% cbba_MRTA() - Solves an MRTA problem with the consensus-based bundle algorithm (CBBA) of Section VII-D
%
% [z, V, info] = cbba_MRTA(prob, a_loc, t_loc, Name=Value)
%
% Writes the problem to a text file (export_GA_instance), solves it with the CBBA executable
% bin/cbba_mrta.exe, reads its solution (read_CBBA_solution) and checks it against the central MILP (9)
% (evaluate_GA_solution). CBBA solves the same problem as PReDuaL and the GA: same cost (8), each pending task
% is performed once by a robot that can perform it, and the route of each robot is a simple path (built from
% its tasks with the decoder of the GA).
%
% INPUTS:
% @prob: Instance of MRTA (with the same costs as for the other solvers, e.g., after prob.normalize_costs())
% @a_loc: Cell array with the location (name) of each agent
% @t_loc: Cell array with the locations (names) of the pending tasks of each task type
%
% NAME/VALUE INPUTS:
% 'exe': Path of the CBBA executable (default: bin/cbba_mrta.exe of this package)
% 'network': Communication network of the robots: "full" (default), "ring", "line", "star" or "random"
% 'prob': Probability of each link of the "random" network (default 0.3; the network is then connected)
% 'seed': Seed of the "random" network (default 1)
% 'bundle': Maximum number of tasks in the bundle of a robot (default 0: no limit)
% 'warp': Bid warping, which guarantees convergence (default true; false: original CBBA)
% 'reorder': After convergence, each robot improves the order of its own path (default false)
% 'maxrounds': Maximum number of rounds (default 10000)
% 'time': Time limit [s] (default inf)
% 'args': Additional options of cbba_mrta (text)
% 'folder': Folder for the instance and solution files (default: tempdir)
% 'name': Name of the instance (default "instance")
% 'keep_files': Keep the instance and solution files (default false)
% 'verbose': Print the output of CBBA (default false)
%
% OUTPUTS:
% @z: Solution of the central MILP (9) given by CBBA ([] if it is not valid)
% @V: Cost (8) of z (Inf if it is not feasible)
% @info: Structure with fields
%   .e_flag: 1 if CBBA found a feasible solution of (9), -1 otherwise
%   .feasible: true if z is a feasible solution of (9)
%   .converged: true if CBBA converged
%   .time: Serial computation time of CBBA [s] (all the robots one after the other; without files)
%   .time_par: Distributed computation time [s]: sum over the rounds of the largest computation time of a
%              robot, as the computation time of PReDuaL (communication delays not included)
%   .time_best, .time_par_best: Serial and distributed time [s] at which the final solution was known: end of
%              the last round that changed the bundles or winners (convergence is detected one round later)
%   .rounds_best: Last round that changed the bundles or winners
%   .wall_time: Wall-clock time of the call to the executable [s] (including starting it and the files)
%   .rounds, .messages: Communication rounds and messages
%   .network, .diameter: Communication network and its diameter
%   .cost_cbba: Cost (8) computed by CBBA (equal to V up to rounding errors)
%   .routes, .tasks: Route (locations) and tasks of each robot (see read_CBBA_solution)
%   .task_v, .task_t: Location and type of each task (see export_GA_instance)
%   .viol, .msg: See evaluate_GA_solution
%   .files: Instance and solution files (deleted unless 'keep_files' is true)
%

function [z, V, info] = cbba_MRTA(prob, a_loc, t_loc, opt)
    arguments
        prob MRTA
        a_loc
        t_loc
        opt.exe {mustBeText} = ""
        opt.network {mustBeMember(opt.network, ["full", "ring", "line", "star", "random"])} = "full"
        opt.prob (1, 1) double {mustBeInRange(opt.prob, 0, 1)} = 0.3
        opt.seed (1, 1) {mustBeNonnegative, mustBeInteger} = 1
        opt.bundle (1, 1) {mustBeNonnegative, mustBeInteger} = 0
        opt.warp (1, 1) {mustBeNumericOrLogical} = true
        opt.reorder (1, 1) {mustBeNumericOrLogical} = false
        opt.maxrounds (1, 1) {mustBePositive, mustBeInteger} = 10000
        opt.time (1, 1) double {mustBePositive} = inf
        opt.args {mustBeText} = ""
        opt.folder {mustBeText} = ""
        opt.name {mustBeText} = "instance"
        opt.keep_files (1, 1) {mustBeNumericOrLogical} = false
        opt.verbose (1, 1) {mustBeNumericOrLogical} = false
    end

    exe = cbba_executable(opt.exe);
    folder = char(opt.folder);
    if isempty(folder); folder = tempdir; end
    if ~isfolder(folder); mkdir(folder); end
    base = tempname(folder); % Unique name
    inst_file = [base, '_instance.txt'];
    sol_file = [base, '_solution.txt'];
    info.files = {inst_file, sol_file};

    % Write the instance (same file as for the GA)
    [task_v, task_t] = export_GA_instance(prob, a_loc, t_loc, inst_file, opt.name);

    % Run CBBA
    args = {inst_file, sol_file, '--network', char(opt.network), '--prob', sprintf('%.17g', opt.prob), ...
            '--seed', sprintf('%d', opt.seed), '--bundle', sprintf('%d', opt.bundle), ...
            '--maxrounds', sprintf('%d', opt.maxrounds)};
    if isfinite(opt.time); args = [args, {'--time', sprintf('%.17g', opt.time)}]; end
    if ~opt.warp; args = [args, {'--no-warp'}]; end
    if opt.reorder; args = [args, {'--reorder'}]; end
    if ~opt.verbose; args = [args, {'--quiet'}]; end
    t_wall = tic;
    [status, out] = system(command_line(exe, args, strtrim(char(opt.args))));
    info.wall_time = toc(t_wall);
    if opt.verbose; fprintf("%s", out); end
    if status ~= 0 || ~isfile(sol_file)
        cleanup_files(info.files, opt.keep_files);
        error("cbba_MRTA:cbba", "CBBA (%s) failed with exit code %d:\n%s", exe, status, out);
    end

    % Read and check the solution
    sol = read_CBBA_solution(sol_file);
    ev = evaluate_GA_solution(prob, a_loc, t_loc, task_v, task_t, sol);
    cleanup_files(info.files, opt.keep_files);

    z = ev.z;
    V = ev.V;
    info.feasible = ev.feasible;
    info.e_flag = 1;
    if ~ev.feasible; info.e_flag = -1; end
    info.converged = sol.converged;
    info.time = sol.time;
    info.time_par = sol.time_par;
    info.time_best = sol.time_best;
    info.time_par_best = sol.time_par_best;
    info.rounds_best = sol.rounds_best;
    info.rounds = sol.rounds;
    info.messages = sol.messages;
    info.network = sol.network;
    info.diameter = sol.diameter;
    info.cost_cbba = sol.cost;
    info.feasible_cbba = sol.feasible;
    info.routes = sol.routes;
    info.tasks = sol.tasks;
    info.task_v = task_v;
    info.task_t = task_t;
    info.viol = ev.viol;
    info.msg = ev.msg;

    if sol.feasible && ~ev.feasible
        warning("cbba_MRTA:check", "The solution of CBBA is not a feasible solution of (9): %s", ev.msg);
    elseif ev.feasible && abs(V - sol.cost) > 1e-6*max(1, abs(V))
        warning("cbba_MRTA:cost", "The cost computed by CBBA (%.10g) differs from the cost (8) of its solution (%.10g)", sol.cost, V);
    end

end

function exe = cbba_executable(exe)
    % Path of the CBBA executable
    exe = char(exe);
    if ~isempty(exe)
        if ~isfile(exe)
            error("cbba_MRTA:exe", "The CBBA executable %s was not found", exe);
        end
        return
    end
    bin = fullfile(fileparts(mfilename('fullpath')), 'bin');
    if ispc
        exe = fullfile(bin, 'cbba_mrta.exe');
    else
        exe = fullfile(bin, 'cbba_mrta'); % Not included in the package (only the Windows executable is)
    end
    if ~isfile(exe)
        error("cbba_MRTA:exe", "The CBBA executable %s was not found (use the option 'exe' to give its path)", exe);
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
