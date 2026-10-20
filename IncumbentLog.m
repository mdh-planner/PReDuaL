%% IncumbentLog - Records the new best integer solutions found by intlinprog (see MIQPholder.solve, track_best)
%
% Used as output function of intlinprog: options = optimoptions(options, 'OutputFcn', @(x, v, s) log.record(x, v, s))
%

classdef IncumbentLog < handle

    properties
        hist = zeros(0, 2) % [time, objective] of each new best integer solution
    end

    properties (Access = private)
        t0 = [] % Start of the solver (if optimValues has no time)
    end

    methods
        function stop = record(self, x, values, state)
            % Output function of intlinprog: x is nonempty when a new integer feasible solution was found
            stop = false;
            if strcmp(state, 'init') || isempty(self.t0); self.t0 = tic; end
            if ~strcmp(state, 'iter') || isempty(x) || ~isfinite(values.fval); return; end
            if isfield(values, 'time') && ~isempty(values.time)
                t = values.time;
            else
                t = toc(self.t0);
            end
            if isempty(self.hist) || values.fval < self.hist(end, 2)
                self.hist(end+1, :) = [t, values.fval];
            end
        end
    end

end
