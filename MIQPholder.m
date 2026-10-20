%% A class for defining, solving and working with MIQP problems
%
% The primal MIQP optimization problem is given by
%
% min_z 0.5*z'*Q*z + c'*z
%  s.t. lb <= A*z <= ub
%  lbz <= z <= ubz
%  C*z = b
%  D_i*z \in {li_i, ui_i}, i = 1,...,p
%
% where dim(z) = n, dim(A) = [m, n], dim(C) = [q, n], dim(D) = [p, n].
%
% Only MILP problems (Q = 0) can be solved, using MOSEK if it is installed and
% MATLAB's intlinprog (Optimization Toolbox) otherwise. See MIQPholder.solve().
%

classdef MIQPholder < handle

    properties
        primal % For storing the primal ingredients
        n % Dimension of primal decision variables
        m % Dimension of primal inequality constraints
        q % Dimension of primal equality constraints
        p % Dimension of number of integer constraints
        nd % Dimension of dual variables
    end

    methods

        % Constructor
        function obj = MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui, lbz, ubz, idx_int)
            % MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui, lbz, ubz, idx_int)
            % The constructor simply takes the ingredients of the MIQP problem
            arguments
                Q (:, :)
                c (:, 1)
                A (:, :)
                lb (:, 1) = []
                ub (:, 1) = []
                C (:, :) = []
                b (:, 1) = []
                D (:, :) = []
                li (:, 1) = []
                ui (:, 1) = []
                lbz (:, 1) = []
                ubz (:, 1) = []
                idx_int (:, 1) = []
            end

            % Get sizes
            obj.n = length(c);
            obj.m = size(A, 1);
            obj.q = size(C, 1);
            obj.p = size(D, 1);
            obj.nd = 2*obj.m + 2*obj.p + obj.q;

            % Give values to empty vectors
            if isempty(Q); Q = sparse(obj.n, obj.n); end
            if isempty(A); A = sparse(0, obj.n); end
            if isempty(C); C = sparse(0, obj.n); end
            if isempty(lb); lb = -inf*ones(obj.m, 1); end
            if isempty(ub); ub = inf*ones(obj.m, 1); end
            if isempty(b); b = zeros(obj.q, 1); end
            if isempty(li); li = zeros(obj.p, 1); end
            if isempty(ui); ui = ones(obj.p, 1); end
            if isempty(lbz); lbz = -inf*ones(obj.n, 1); end
            if isempty(ubz); ubz = inf*ones(obj.n, 1); end
            if isempty(idx_int) && ~isempty(D)
                idx_int = find(any(D, 1))';
            end

            % Assign variables
            obj.primal.Q = Q;
            obj.primal.c = c;
            obj.primal.A = A;
            obj.primal.lb = lb;
            obj.primal.ub = ub;
            obj.primal.C = C;
            obj.primal.b = b;
            obj.primal.D = D;
            obj.primal.li = li;
            obj.primal.ui = ui;
            obj.primal.lbz = lbz;
            obj.primal.ubz = ubz;
            obj.primal.idx_int = idx_int;

        end

        function [z_out, V_out, info] = solve(self, opt)
            % [z, V, info] = MIQPholder.solve(Name=Value) - Solves the MILP problem
            %
            % NAME/VALUE INPUTS:
            %   solver: "auto" (MOSEK if installed, else intlinprog), "mosek" or "intlinprog"
            %   max_time: Maximum solver time in seconds (default inf)
            %   verbose: Shows solver output (default false)
            %   mosek_param: Structure with additional MOSEK parameters (MSK_*)
            %   warmstart: Initial integer solution (used by MOSEK, and by intlinprog from MATLAB R2020b)
            %   relax: Solve the LP relaxation (integrality constraints are ignored) and return its
            %          dual variables. Uses MOSEK or linprog (default false)
            %   track_best: Record when the solver found each new best integer solution (default false).
            %               For MOSEK, the log of its mixed-integer optimizer (written to a temporary file and
            %               captured from the screen output) is parsed: its lines are printed when new solutions
            %               are found and every MSK_IPAR_LOG_MIO_FREQ relaxations, and its TIME column has a
            %               resolution of 0.1 s or 0.01 s. For intlinprog, an output function is used
            %
            % OUTPUTS:
            %   @z_out: Solution (NaN if no solution was found)
            %   @V_out: Objective value c'*z_out
            %   @info: Structure with fields
            %       .e_flag: 1 optimal, 2 feasible but not proven optimal (e.g. time limit), -1 no solution
            %       .feasible: true if a solution was returned
            %       .obj_bound: Lower bound of the optimal value. If e_flag == 1, the bound provided by the
            %                   solver (or V_out minus the default gap tolerance if it is not available);
            %                   otherwise, the bound provided by the solver when it stopped, e.g., at the
            %                   time limit (-Inf if not available)
            %       .y_eq, .y_in: (only if relax = true) Dual variables of the constraints C*z = b and
            %                     lb <= A*z <= ub, such that the reduced costs are c - C'*y_eq - A'*y_in
            %       .run_time: Solver time [s] (MOSEK optimizer time, or wall time for intlinprog/linprog)
            %       .wall_time: Wall time of the call [s]
            %       .solver: Name of the solver that was used
            %       .time_best: (track_best) Solver time [s] at which the returned solution was found, on the
            %                   clock of .run_time (NaN if it could not be determined)
            %       .best_hist: (track_best) [time, objective] of each new best integer solution
            arguments
                self MIQPholder
                opt.solver (1, 1) string {mustBeMember(opt.solver, ["auto", "mosek", "intlinprog"])} = "auto"
                opt.max_time (1, 1) double {mustBePositive} = inf
                opt.verbose (1, 1) {mustBeNumericOrLogical} = false
                opt.mosek_param struct = struct()
                opt.warmstart (:, 1) double = []
                opt.relax (1, 1) {mustBeNumericOrLogical} = false
                opt.track_best (1, 1) {mustBeNumericOrLogical} = false
            end
            opt.track_best = opt.track_best && ~opt.relax;

            if nnz(self.primal.Q) > 0
                error("MIQPholder:not_MILP", "MIQPholder.solve() can only solve MILP problems (Q = 0)");
            end

            solver = opt.solver;
            if strcmp(solver, "auto")
                solver = MIQPholder.default_solver();
            end

            start_time = tic;
            switch solver
                case "mosek"
                    [z_out, info] = self.solve_with_mosek(opt);
                case "intlinprog"
                    [z_out, info] = self.solve_with_intlinprog(opt);
            end

            % If MOSEK fails for a reason other than infeasibility, use intlinprog/linprog instead
            % (only if the solver was selected automatically)
            if strcmp(solver, "mosek") && info.e_flag < 0 && strcmp(opt.solver, "auto") && MIQPholder.has_intlinprog() ...
                    && ~contains(info.mosek.prosta, 'INFEAS')
                MIQPholder.warn_fallback(info);
                info_mosek = info;
                [z_out, info] = self.solve_with_intlinprog(opt);
                info.mosek_failure = info_mosek;
                solver = "intlinprog";
            end
            info.wall_time = toc(start_time);
            if ~isfield(info, 'run_time') || isnan(info.run_time)
                info.run_time = info.wall_time;
            end
            info.solver = solver;
            info.feasible = info.e_flag > 0;

            if info.feasible
                V_out = self.primal.c'*z_out;
            else
                z_out = NaN(self.n, 1);
                V_out = NaN;
            end

            % Lower bound of the optimal value. For an optimal solution, use the bound of the solver if it
            % is consistent, and otherwise the default relative gap tolerance of the solvers (1e-4). If the
            % solver stopped before proving optimality, keep its bound (-Inf if none). For LPs, the bound
            % is the optimal value.
            bound = NaN;
            if isfield(info, 'obj_bound'); bound = info.obj_bound; end
            if opt.relax
                bound = V_out;
            elseif info.e_flag == 1
                if ~(isfinite(bound) && bound <= V_out + 1e-9 && bound >= V_out - 1e-3*max(1, abs(V_out)))
                    bound = V_out - 1e-4*max(1, abs(V_out));
                end
            elseif ~(isfinite(bound) && ~(bound > V_out + 1e-9)) % V_out is NaN if there is no solution
                bound = -inf;
            end
            info.obj_bound = bound;

            % Time at which the returned solution was found: first new best solution with its cost
            if opt.track_best
                if ~isfield(info, 'best_hist'); info.best_hist = []; end
                info.time_best = MIQPholder.time_of_best(info.best_hist, V_out, info.run_time);
            end

        end

        function to_sparse(self)
            % Converts the main ingredients of the MIQP problem to sparse matrices
            self.primal.Q = sparse(self.primal.Q);
            self.primal.A = sparse(self.primal.A);
            self.primal.C = sparse(self.primal.C);
            self.primal.D = sparse(self.primal.D);
        end

        function [feas, viol] = is_feasible(self, z, check_int, tol)
            % [feas, viol] = is_feasible(z, [check_int], [tol])
            % Returns true if @z satisfies all the constraints of the problem (and the integrality
            % constraints if @check_int is true), with tolerance @tol (relative to the right-hand sides).
            % @viol is the largest (scaled) constraint violation.
            arguments
                self MIQPholder
                z (:, 1) double
                check_int (1, 1) {mustBeNumericOrLogical} = true
                tol (1, 1) double = 1e-5
            end
            p = self.primal;
            if length(z) ~= self.n || any(~isfinite(z))
                feas = false;
                viol = inf;
                return
            end
            scaled = @(v, rhs) max([0; v(:)./(1 + abs(rhs(:)))]); % Largest violation scaled by the right-hand side
            Az = p.A*z;
            viol = max([scaled(abs(p.C*z - p.b), p.b), scaled(p.lb - Az, p.lb), scaled(Az - p.ub, p.ub), ...
                        scaled(p.lbz - z, p.lbz), scaled(z - p.ubz, p.ubz)]);
            feas = viol <= tol;
            if check_int && ~isempty(p.idx_int)
                viol_int = max([0; abs(z(p.idx_int) - round(z(p.idx_int)))]);
                feas = feas && viol_int <= 1e-4;
                viol = max(viol, viol_int);
            end
        end

    end % End of public methods

    methods (Access = private)

        function [z_out, info] = solve_with_mosek(self, opt)
            % Solves the MILP problem using MOSEK's mosekopt
            if ~MIQPholder.has_mosek()
                error("MIQPholder:no_mosek", "MOSEK (mosekopt) was not found on the MATLAB path");
            end

            prob.c = full(self.primal.c);
            prob.a = sparse([self.primal.C; self.primal.A]);
            prob.blc = [self.primal.b; self.primal.lb];
            prob.buc = [self.primal.b; self.primal.ub];
            prob.blx = self.primal.lbz;
            prob.bux = self.primal.ubz;
            if ~opt.relax
                prob.ints.sub = self.primal.idx_int(:)';
            end
            if ~isempty(opt.warmstart) && ~opt.relax
                prob.sol.int.xx = opt.warmstart;
            end

            param = opt.mosek_param;
            if isfinite(opt.max_time)
                param.MSK_DPAR_OPTIMIZER_MAX_TIME = opt.max_time;
            end
            if ~isfield(param, 'MSK_IPAR_LOG')
                param.MSK_IPAR_LOG = double(opt.verbose);
            end
            cmd = 'minimize info';
            if opt.track_best
                % The log of the mixed-integer optimizer has the time of each new best solution. It is written
                % to a temporary file and also captured from the screen output, and the first one that has it
                % is parsed
                param.MSK_IPAR_LOG = max(double(opt.verbose), 10);
                param.MSK_IPAR_LOG_MIO = 4;
                log_file = [tempname, '.log'];
                if any(isspace(log_file))
                    log_file = '';
                else
                    cmd = [cmd, ' log(', log_file, ')'];
                end
                r = -1; res = struct(); %#ok<NASGU> Assigned by evalc
                screen = evalc('[r, res] = mosekopt(cmd, prob, param);');
                if opt.verbose; fprintf('%s', screen); end
                info.best_hist = [];
                if ~isempty(log_file) && isfile(log_file)
                    info.best_hist = MIQPholder.parse_mosek_mio_log(fileread(log_file));
                    delete(log_file);
                end
                if isempty(info.best_hist)
                    info.best_hist = MIQPholder.parse_mosek_mio_log(screen);
                end
            else
                if ~opt.verbose; cmd = [cmd, ' echo(0)']; end
                [r, res] = mosekopt(cmd, prob, param);
            end

            info.mosek_e_flag = r;
            info.mosek.rcodestr = '';
            if isfield(res, 'rcodestr'); info.mosek.rcodestr = res.rcodestr; end
            info.mosek.prosta = 'UNKNOWN';
            info.mosek.solsta = 'UNKNOWN';
            z_out = NaN(self.n, 1);
            sol = struct();
            if isfield(res, 'sol')
                if opt.relax && isfield(res.sol, 'bas')
                    sol = res.sol.bas; % Basic solution (LP)
                elseif opt.relax && isfield(res.sol, 'itr')
                    sol = res.sol.itr; % Interior-point solution (LP)
                elseif ~opt.relax && isfield(res.sol, 'int')
                    sol = res.sol.int; % Integer solution (MILP)
                end
            end
            if isfield(sol, 'prosta'); info.mosek.prosta = upper(char(string(sol.prosta))); end
            if isfield(sol, 'solsta'); info.mosek.solsta = upper(char(string(sol.solsta))); end
            has_x = isfield(sol, 'xx') && length(sol.xx) == self.n && all(isfinite(sol.xx));
            info.mosek.violation = inf;
            if has_x
                info.mosek.sol = sol;
                z_out = sol.xx(:);
                [sol_feasible, info.mosek.violation] = self.is_feasible(z_out, ~opt.relax);
            end

            % Accept the solution if it satisfies the constraints. It is optimal if MOSEK terminated
            % normally (r = 0) or reports an optimal solution status.
            if has_x && sol_feasible
                if r == 0 || contains(info.mosek.solsta, 'OPTIMAL')
                    info.e_flag = 1; % Optimal solution found
                else
                    info.e_flag = 2; % Feasible solution (e.g., time limit reached)
                end
            else
                info.e_flag = -1; % No solution (infeasible, out of memory, time limit without solution, etc.)
            end

            if opt.relax && isfield(sol, 'y')
                q = size(self.primal.C, 1);
                info.y_eq = sol.y(1:q);
                info.y_in = sol.y(q+1:end);
            end

            info.run_time = NaN;
            info.obj_bound = NaN;
            if isfield(res, 'info')
                info.mosek.info = res.info;
                if isfield(res.info, 'MSK_DINF_OPTIMIZER_TIME')
                    info.run_time = res.info.MSK_DINF_OPTIMIZER_TIME;
                end
                if ~opt.relax && isfield(res.info, 'MSK_DINF_MIO_OBJ_BOUND')
                    info.obj_bound = res.info.MSK_DINF_MIO_OBJ_BOUND;
                end
            end
        end

        function [z_out, info] = solve_with_intlinprog(self, opt)
            % Solves the MILP problem using MATLAB's intlinprog (Optimization Toolbox)
            if ~MIQPholder.has_intlinprog()
                error("MIQPholder:no_intlinprog", "intlinprog (Optimization Toolbox) was not found");
            end

            % Write lb <= A*z <= ub as Ain*z <= bin
            A = sparse(self.primal.A);
            fin_ub = isfinite(self.primal.ub);
            fin_lb = isfinite(self.primal.lb);
            Ain = [A(fin_ub, :); -A(fin_lb, :)];
            bin = [self.primal.ub(fin_ub); -self.primal.lb(fin_lb)];

            if opt.verbose; disp_level = 'iter'; else; disp_level = 'off'; end

            if opt.relax
                % LP relaxation using linprog
                options = MIQPholder.optim_options('linprog', disp_level, inf);
                t_start = tic;
                [z_out, ~, exitflag, output, lambda] = MIQPholder.call_matlab_optim('linprog', full(self.primal.c), Ain, bin, ...
                    sparse(self.primal.C), self.primal.b, self.primal.lbz, self.primal.ubz, options);
                info.run_time = toc(t_start);
                info.linprog.exitflag = exitflag;
                info.linprog.output = output;
                if exitflag == 1
                    info.e_flag = 1;
                    % linprog's multipliers satisfy c + C'*eqlin + Ain'*ineqlin - lower + upper = 0
                    n_ub = nnz(fin_ub);
                    info.y_eq = -lambda.eqlin;
                    info.y_in = zeros(size(self.primal.A, 1), 1);
                    info.y_in(fin_ub) = -lambda.ineqlin(1:n_ub);
                    info.y_in(fin_lb) = info.y_in(fin_lb) + lambda.ineqlin(n_ub+1:end);
                else
                    info.e_flag = -1;
                end
                if isempty(z_out); z_out = NaN(self.n, 1); end
                return
            end

            max_time = opt.max_time;
            if ~isfinite(max_time); max_time = 1e7; end
            options = MIQPholder.optim_options('intlinprog', disp_level, max_time);
            if opt.track_best
                incumbents = IncumbentLog();
                options = optimoptions(options, 'OutputFcn', @(x, values, state) incumbents.record(x, values, state));
            end
            args = {full(self.primal.c), self.primal.idx_int(:)', Ain, bin, sparse(self.primal.C), self.primal.b, ...
                    self.primal.lbz, self.primal.ubz};
            if ~isempty(opt.warmstart) && length(opt.warmstart) == self.n && MIQPholder.intlinprog_has_x0()
                args{end+1} = opt.warmstart; % Initial point x0
            end

            t_start = tic;
            [z_out, ~, exitflag, output] = MIQPholder.call_matlab_optim('intlinprog', args{:}, options);
            info.run_time = toc(t_start);
            info.intlinprog.exitflag = exitflag;
            info.intlinprog.output = output;
            if opt.track_best
                info.best_hist = incumbents.hist;
            end

            info.obj_bound = NaN;
            if isstruct(output) && isfield(output, 'absolutegap') && ~isempty(output.absolutegap) && ~isempty(z_out)
                info.obj_bound = full(self.primal.c)'*z_out - output.absolutegap;
            end
            if exitflag == 1 || exitflag == 3
                info.e_flag = 1; % Optimal solution found
            elseif exitflag == 2
                info.e_flag = 2; % Feasible solution, stopped prematurely (e.g., time limit)
            else
                info.e_flag = -1; % No solution
            end
            if isempty(z_out); z_out = NaN(self.n, 1); end
        end

    end % End of private methods

    methods (Static)

        function solver = default_solver()
            % Returns the solver used by solver="auto": MOSEK if installed, otherwise intlinprog
            if MIQPholder.has_mosek()
                solver = "mosek";
            elseif MIQPholder.has_intlinprog()
                solver = "intlinprog";
            else
                error("MIQPholder:no_solver", ...
                    "No MILP solver found. Install MOSEK (https://docs.mosek.com) or the Optimization Toolbox (intlinprog).");
            end
        end

        function warn_fallback(info)
            % Warns (once per MATLAB session) that MOSEK failed and intlinprog is used instead
            persistent warned
            if isempty(warned)
                warned = true;
                warning("MIQPholder:mosek_failed", ...
                    "MOSEK could not solve a problem (response code %d: %s, problem status: %s, solution status: %s, " + ...
                    "max. constraint violation of the returned solution: %g). Using intlinprog/linprog for the " + ...
                    "problems where MOSEK fails. This warning is shown only once.", ...
                    info.mosek_e_flag, info.mosek.rcodestr, info.mosek.prosta, info.mosek.solsta, info.mosek.violation);
            end
        end

        function flag = has_mosek()
            % True if MOSEK is on the path and has a valid license (checked once per MATLAB session
            % by solving a tiny LP; run "clear MIQPholder" after renewing the license)
            persistent ok
            if isempty(ok)
                ok = exist('mosekopt', 'file') > 0;
                if ok
                    prob = struct('c', 1, 'a', sparse(1), 'blc', 0, 'buc', 1, 'blx', 0, 'bux', 1);
                    try
                        [r, res] = mosekopt('minimize echo(0)', prob);
                        ok = r == 0;
                        msg = '';
                        if isfield(res, 'rcodestr'); msg = res.rcodestr; end
                    catch ME
                        ok = false;
                        msg = ME.message;
                    end
                    if ~ok
                        warning("MIQPholder:mosek_unavailable", ...
                            "MOSEK is installed but cannot be used (%s). Using intlinprog/linprog instead.", msg);
                    end
                end
            end
            flag = ok;
        end

        function varargout = call_matlab_optim(fname, varargin)
            % Calls the Optimization Toolbox function @fname (intlinprog or linprog). The MOSEK toolbox
            % ships functions with the same names, which are temporarily removed from the path (unless
            % optim_path_guard already did it, since changing the path is slow).
            p = which(fname);
            if ~isempty(p) && ~startsWith(p, matlabroot)
                p0 = path;
                restore = onCleanup(@() path(p0)); %#ok<NASGU>
                p = MIQPholder.remove_shadowing({fname});
                p = p{1};
            end
            if isempty(p)
                error("MIQPholder:no_optim", "%s (Optimization Toolbox) was not found", fname);
            end
            [varargout{1:nargout}] = feval(fname, varargin{:});
        end

        function guard = optim_path_guard(solver_args)
            % guard = optim_path_guard(solver_args)
            % If the problems solved with the options @solver_args of MIQPholder.solve use intlinprog/linprog
            % and other functions with the same names shadow them (e.g., those of the MOSEK toolbox), their
            % folders are removed from the path until the returned object is destroyed. Returns [] otherwise.
            guard = [];
            solver = "auto";
            for q = 1:2:(length(solver_args) - 1)
                if strcmpi(string(solver_args{q}), "solver"); solver = string(solver_args{q+1}); end
            end
            mosek_ok = MIQPholder.has_mosek(); % Check (and cache) it before changing the path
            if strcmp(solver, "mosek") || (strcmp(solver, "auto") && mosek_ok)
                return
            end
            p0 = path;
            [~, changed] = MIQPholder.remove_shadowing({'intlinprog', 'linprog'});
            if changed
                guard = onCleanup(@() path(p0));
            end
        end

        function [p, changed] = remove_shadowing(fnames)
            % Removes from the path the folders whose functions @fnames shadow those of MATLAB (not under
            % matlabroot). Returns the resolved file of each function and true if the path was changed.
            changed = false;
            p = cell(size(fnames));
            for q = 1:numel(fnames)
                p{q} = which(fnames{q});
                while ~isempty(p{q}) && ~startsWith(p{q}, matlabroot)
                    rmpath(fileparts(p{q}));
                    changed = true;
                    p{q} = which(fnames{q});
                end
            end
        end

        function options = optim_options(fname, disp_level, max_time)
            % Options of intlinprog/linprog (cached, since optimoptions is slow compared to small MILPs)
            persistent cache
            if isempty(cache); cache = struct('key', {}, 'options', {}); end
            key = sprintf('%s|%s|%.17g', fname, disp_level, max_time);
            q = find(strcmp({cache.key}, key), 1);
            if isempty(q)
                if strcmp(fname, 'intlinprog')
                    options = optimoptions('intlinprog', 'Display', disp_level, 'MaxTime', max_time);
                else
                    options = optimoptions('linprog', 'Display', disp_level);
                end
                cache(end+1) = struct('key', key, 'options', options);
            else
                options = cache(q).options;
            end
        end

        function hist = parse_mosek_mio_log(txt)
            % Returns the rows [time, objective] of the log @txt of MOSEK's mixed-integer optimizer in which the
            % best integer solution is known (columns TIME and BEST_INT_OBJ of its table), in the order of the
            % log, or [] if the log has no such table (e.g., the log was not written)
            hist = [];
            col_obj = 0;
            for line = regexp(txt, '\r?\n', 'split')
                tok = strsplit(strtrim(line{1}));
                if any(strcmp(tok, 'BEST_INT_OBJ')) && any(strcmp(tok, 'TIME')) % Header (repeated after a restart)
                    col_obj = find(strcmp(tok, 'BEST_INT_OBJ'), 1);
                    col_time = find(strcmp(tok, 'TIME'), 1);
                    num_cols = numel(tok);
                    if isempty(hist); hist = zeros(0, 2); end
                    continue
                end
                if col_obj == 0 || numel(tok) ~= num_cols; continue; end
                vals = str2double(tok); % "NA" (no integer solution yet) gives NaN
                if isnan(vals(1)) || isnan(vals(col_time)); continue; end % Not a row of the table
                if ~isnan(vals(col_obj))
                    hist(end+1, :) = [vals(col_time), vals(col_obj)]; %#ok<AGROW>
                end
            end
        end

        function t = time_of_best(hist, V, t_end)
            % Time of the first row of @hist ([time, objective] of the new best solutions) whose objective is
            % at most V (the cost of the returned solution, up to the precision of the log). If no row matches
            % (the solution was found after the last recorded one), the solver time @t_end. NaN if @V is not
            % finite, or if @hist is [] (nothing recorded, e.g., the log could not be read)
            t = NaN;
            if ~isfinite(V); return; end
            if isequal(size(hist), [0, 0])
                warning("MIQPholder:time_best", ['The time of the best solution could not be determined ', ...
                    '(the solver log has no table of integer solutions)']);
                return
            end
            q = find(hist(:, 2) <= V + 1e-6*max(1, abs(V)), 1);
            if isempty(q)
                t = t_end;
            else
                t = hist(q, 1);
            end
        end

        function flag = intlinprog_has_x0()
            % True if intlinprog accepts an initial point x0 (MATLAB R2020b or later)
            persistent ok
            if isempty(ok)
                ok = exist('isMATLABReleaseOlderThan', 'file') > 0; % Introduced in R2020b
            end
            flag = ok;
        end

        function flag = has_intlinprog()
            % True if the Optimization Toolbox versions of intlinprog and linprog are available (MOSEK ships
            % functions with the same names). Checked once per session (run "clear MIQPholder" to reset).
            persistent ok
            if isempty(ok)
                in_toolbox = @(f) any(startsWith(cellstr(which(f, '-all')), matlabroot));
                ok = license('test', 'Optimization_Toolbox') && in_toolbox('intlinprog') && in_toolbox('linprog');
            end
            flag = ok;
        end

    end % End of static methods

end
