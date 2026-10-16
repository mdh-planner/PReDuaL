%% A class for solving the MILP problems of the agents, sequentially or in parallel
%
% Solves batches of MILP problems of the agents: the problem of agent i has the constraints of its
% MIQPholder (set Z_i in (11)) and the given costs c and bounds lbz, ubz (the default bounds of the
% agent if not given). It is used for the problems (13) of Algorithm 1 and for those of the primal
% recovery.
%
% If 'parallel' is true, the problems of each batch are solved in parallel with parfor (Parallel
% Computing Toolbox), as in the distributed setting of the article, where each robot solves its own
% problem. The constraint matrices are sent to the workers once (parallel.pool.Constant), and then
% only the costs and bounds are sent. If the parallel pool cannot be used, the problems are solved
% sequentially.
%

classdef AgentSolver < handle

    properties
        miqp % 1 x m cell with the MIQPholder of each agent
        lbz % 1 x m cell with the default lower bounds of the decision variables of each agent
        ubz % 1 x m cell with the default upper bounds of the decision variables of each agent
        solver_args % Cell array with name/value options for MIQPholder.solve
        lp_first = false % If true, the LP relaxation is solved first (see solve_one)
        parallel = false % true if the problems are solved in parallel
        const % parallel.pool.Constant with copies of the MIQPholder objects (used by the workers)
    end

    methods

        function obj = AgentSolver(mrta, solver_args, use_parallel)
            % AgentSolver(mrta, [solver_args], [use_parallel])
            %   @mrta: Instance of MRTA with the agents built by build_distMILP
            %   @solver_args: Cell array with name/value options for MIQPholder.solve
            %   @use_parallel: Solve the problems in parallel (starts a parallel pool if needed)
            arguments
                mrta MRTA
                solver_args cell = {}
                use_parallel (1, 1) {mustBeNumericOrLogical} = false
            end
            m = mrta.m;
            obj.miqp = cell(1, m);
            obj.lbz = cell(1, m);
            obj.ubz = cell(1, m);
            for i = 1:m
                obj.miqp{i} = mrta.agent{i}.miqp;
                obj.lbz{i} = obj.miqp{i}.primal.lbz;
                obj.ubz{i} = obj.miqp{i}.primal.ubz;
            end
            obj.solver_args = solver_args;
            if use_parallel
                obj.parallel = obj.start_parallel();
            end
        end

        function [z, V, info, wall_time] = solve(self, agents, c, lbz, ubz, ws)
            % [z, V, info, wall_time] = solve(agents, c, [lbz], [ubz], [ws])
            % Solves the MILP problems of the agents @agents (vector of indices, an agent can appear more
            % than once) with the costs c{q} and bounds lbz{q}, ubz{q} (default bounds if empty).
            % ws{q} is a feasible solution used to warm-start the solver (none if empty).
            % Returns cell arrays with the solution @z and the information @info of each problem
            % (see MIQPholder.solve), the objective values @V and the wall time of the batch.
            n = length(agents);
            if nargin < 4 || isempty(lbz); lbz = cell(1, n); end
            if nargin < 5 || isempty(ubz); ubz = cell(1, n); end
            if nargin < 6 || isempty(ws); ws = cell(1, n); end
            for q = 1:n
                if isempty(lbz{q}); lbz{q} = self.lbz{agents(q)}; end
                if isempty(ubz{q}); ubz{q} = self.ubz{agents(q)}; end
            end

            start_time = tic;
            if self.parallel && n > 1
                try
                    [z, V, info] = AgentSolver.solve_parfor(self.const, agents, c, lbz, ubz, ws, self.solver_args, self.lp_first);
                catch ME
                    warning("AgentSolver:parallel_failed", ...
                        "The parallel solution of the MILP problems failed (%s). They are solved sequentially from now on.", ME.message);
                    self.parallel = false;
                    [z, V, info] = self.solve_serial(agents, c, lbz, ubz, ws);
                end
            else
                [z, V, info] = self.solve_serial(agents, c, lbz, ubz, ws);
            end
            wall_time = toc(start_time);
        end

        function [z, V, info] = solve_serial(self, agents, c, lbz, ubz, ws)
            % Solves the problems one after the other (see solve)
            n = length(agents);
            z = cell(1, n);
            V = zeros(1, n);
            info = cell(1, n);
            for q = 1:n
                P = self.miqp{agents(q)};
                primal = P.primal; % Restored after solving
                P.primal.c = c{q};
                P.primal.lbz = lbz{q};
                P.primal.ubz = ubz{q};
                args = [self.solver_args, AgentSolver.ws_args(ws{q})];
                [z{q}, V(q), info{q}] = AgentSolver.solve_one(P, args, self.lp_first);
                P.primal = primal;
            end
        end

    end

    methods (Access = private)

        function ok = start_parallel(self)
            % Prepares the parallel pool. Returns false (and the problems are solved sequentially)
            % if the Parallel Computing Toolbox is not available or the workers cannot use the solver.
            ok = false;
            if ~(exist('parpool', 'file') > 0 && license('test', 'Distrib_Computing_Toolbox'))
                warning("AgentSolver:no_parallel_toolbox", ...
                    "The Parallel Computing Toolbox is not available. The MILP problems are solved sequentially.");
                return
            end
            try
                pool = gcp('nocreate');
                if isempty(pool)
                    pool = parpool();
                end

                % Make this package and MOSEK (if it can be used) available on the workers
                dirs = {fileparts(mfilename('fullpath'))};
                if MIQPholder.has_mosek()
                    dirs{end+1} = fileparts(which('mosekopt'));
                end
                wait(parfevalOnAll(pool, @addpath, 0, dirs{:}));

                % The workers must be able to use the same MILP solver as this MATLAB session
                solver = AgentSolver.requested_solver(self.solver_args);
                can_use = fetchOutputs(parfevalOnAll(pool, @AgentSolver.can_use_solver, 1, solver));
                if ~all(can_use)
                    warning("AgentSolver:workers_solver", ...
                        "The workers of the parallel pool cannot use the MILP solver %s (e.g., license or path). " + ...
                        "The MILP problems are solved sequentially.", solver);
                    return
                end

                % Copies of the agents' problems on the workers (sent only once)
                self.const = parallel.pool.Constant(self.miqp);
                ok = true;
            catch ME
                warning("AgentSolver:parallel_setup", ...
                    "The parallel pool could not be set up (%s). The MILP problems are solved sequentially.", ME.message);
            end
        end

    end

    methods (Static)

        function [z, V, info] = solve_parfor(const, agents, c, lbz, ubz, ws, solver_args, lp_first)
            % Solves the problems in parallel. @const.Value is a cell array with copies of the
            % MIQPholder objects of the agents on each worker.
            n = length(agents);
            z = cell(1, n);
            V = zeros(1, n);
            info = cell(1, n);
            parfor q = 1:n
                P = const.Value{agents(q)};
                P.primal.c = c{q};
                P.primal.lbz = lbz{q};
                P.primal.ubz = ubz{q};
                args = [solver_args, AgentSolver.ws_args(ws{q})];
                [z{q}, V(q), info{q}] = AgentSolver.solve_one(P, args, lp_first);
            end
        end

        function [z, V, info] = solve_one(P, args, lp_first)
            % Solves the MILP problem of the MIQPholder @P with the name/value options @args (see
            % MIQPholder.solve). If @lp_first is true, its LP relaxation is solved first with the simplex
            % method: if the solution is integer, it is optimal and the MILP solver is not called.
            % Otherwise, the MILP problem is solved, and the time of the LP is added to info.run_time.
            t_lp = 0;
            if lp_first
                largs = AgentSolver.lp_args(args);
                [z_lp, ~, info_lp] = P.solve(largs{:}, 'relax', true);
                t_lp = info_lp.run_time;
                if info_lp.feasible
                    idx = P.primal.idx_int;
                    z_r = z_lp;
                    z_r(idx) = round(z_lp(idx));
                    if max([0; abs(z_lp(idx) - z_r(idx))]) <= 1e-6 && P.is_feasible(z_r, true)
                        z = z_r;
                        V = P.primal.c'*z;
                        info = info_lp;
                        info.e_flag = 1;
                        info.obj_bound = V;
                        info.lp_integral = true;
                        return
                    end
                end
            end
            [z, V, info] = P.solve(args{:});
            info.run_time = info.run_time + t_lp;
            info.lp_integral = false;
        end

        function args = lp_args(args)
            % Name/value options of MIQPholder.solve for the LP relaxation of solve_one: no warm start, and
            % the dual simplex method for MOSEK (unless MSK_IPAR_OPTIMIZER is given)
            keep = true(1, length(args));
            mosek_param = struct();
            for q = 1:2:(length(args) - 1)
                name = lower(char(string(args{q})));
                if strcmp(name, 'warmstart') || strcmp(name, 'mosek_param')
                    keep([q, q+1]) = false;
                    if strcmp(name, 'mosek_param'); mosek_param = args{q+1}; end
                end
            end
            if ~isfield(mosek_param, 'MSK_IPAR_OPTIMIZER')
                mosek_param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_DUAL_SIMPLEX';
            end
            args = [args(keep), {'mosek_param', mosek_param}];
        end

        function args = ws_args(ws)
            % Name/value options of MIQPholder.solve for the initial solution @ws (none if empty)
            if isempty(ws)
                args = {};
            else
                args = {'warmstart', ws(:)};
            end
        end

        function solver = requested_solver(solver_args)
            % Returns the MILP solver that MIQPholder.solve uses with the options @solver_args
            solver = "auto";
            for q = 1:2:(length(solver_args) - 1)
                if strcmpi(string(solver_args{q}), "solver")
                    solver = string(solver_args{q+1});
                end
            end
            if strcmp(solver, "auto")
                solver = MIQPholder.default_solver();
            end
        end

        function flag = can_use_solver(solver)
            % Returns true if this MATLAB session (e.g., a worker) can use the MILP solver @solver
            switch string(solver)
                case "mosek"
                    flag = MIQPholder.has_mosek();
                case "intlinprog"
                    flag = MIQPholder.has_intlinprog();
                otherwise
                    flag = false;
            end
        end

    end

end
