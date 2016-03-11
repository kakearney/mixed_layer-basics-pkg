function [t,y, varargout] = odewrap(solver, fun, tspan, y0, options, varargin)
%ODEWRAP Wrapper for ODE solvers with more flexible input and output
%
% [t,y] = odewrap(solver, fun, tspan, y0)
% [t,y] = odewrap(solver, fun, tspan, y0, options)
% [t,y] = odewrap(solver, fun, tspan, y0, options, in1, in2, ...)
% [t,y] = odewrap(solver, fun, tspan, y0, [],      in1, in2, ...)
% [t,y, dy, out1, out2, ...] = odewrap(...)
%
% This function is a wrapper for Matlab's ODE solvers.  It adds the
% flexibility of allowing the solvers to evaluate functions of the form
%
% [dydt, out1, out2, ...] = odefun(t, y, in1, in2, ...), 
%
% where both y and dydt can be matrices rather than vectors.  It can be
% used to run any of the variable-step solvers provided with Matlab, as
% well as the fixed-step solvers that can downloaded from the Mathworks
% website (see "Tech Note 1510: Differential Equations in Matlab").
%
% This function is not intended to support systems that require mass matrix
% or Jacobian properties, or those that utilize Events.
%
% Input variables:
%
%   solver:     function handle to ODE solver
%
%   fun:        function handle that evaluates the differential equation.
%               This function should be of the form dydt = fun(t, y), where
%               t is a scalar time value, and y and dydt are matrices of
%               identical size.  See params input to pass additional input
%               to the differential equation.
%
%   tspan:      vector specifying interval of differentiation.  If a
%               variable-step solver is used and tspan includes two
%               elements [t0 tf], the solver returns the solution evaluated
%               at every integration step.  Otherwise, a solution will be
%               returned at each specified time value.  For fixed-step
%               solvers, values will only be returned at the specified
%               values, regardless of the  length of tspan. 
%
%   y0:         matrix of initial conditions.
%
%   options:    structure of optional parameters that change the default
%               integration properties.  Only applicable to variable-step
%               solvers (use empty array if you need to pass additional
%               parameters to a fixed-step solver).  This function is not
%               designed to solve systems where mass matrix or Jacobian
%               properties are needed.  Those properties may work, but if
%               so it is by accident.
%
%   in#:        additional parameters required by fun.  These can be any
%               size or data type, depending on the specific function being
%               evaluated.  Parameters are held constant throughout the
%               integration time span.
%
% Output variables:
%
%   t:          vector of times values corresponding to solution
%
%   y:          array of solutions.  This will be a length(t) x size(y0)
%               array.
%
%   dy:         array of dy/dt values at each of the solution times
%
%   out#:       additional output variables returned by the differential
%               equation function.  These may be any size or data type,
%               depending of the specific function being evaluated.

% Copyright 2008 Kelly Kearney

%------------------------------
% Check input and output
%------------------------------

if nargin(fun) < 2
    error('The differential equation function must be of the form dydt = fun(t,y)');
end

nvar = nargin(fun) - 2;
if nvar ~= length(varargin)
    error('Number of parameters provided does not match number needed by function');
end

nvarout = nargout(fun);
nextraout = nargout - 2;
if nextraout > nvarout & nvarout ~= -1
    error('You have requested more output variables than your function provides');
end

%------------------------------
% Create function that matches
% Matlab's ODE format (one 
% vector ouput, input of scalar 
% time and vector state 
% variable)
%------------------------------

sz = size(y0);
odefun = @(t,y) newfun(t, y, fun, sz, varargin{:});

%------------------------------
% Integrate
%------------------------------

if nargout(solver) == -1    % Adaptive step solvers
    if nargin < 5 || isempty(options)
        [t,y] = feval(solver, odefun, tspan, y0);
    else
        % TODO Figure out which options need modifications, if any
        [t,y] = feval(solver, odefun, tspan, y0, options);
    end
elseif nargout(solver) == 1 % Fixed step solvers
    y = feval(solver, odefun, tspan, y0);
    t = tspan;
elseif nargout(solver) == 2 % ode4nonneg and possibly others
    [t,y] = feval(solver, odefun, tspan, y0);
end
    
%------------------------------
% Reshape results to match 
% input
%------------------------------

nt = length(t);
y = reshape(y, [nt sz]);

%------------------------------
% Rerun ODE to get additional
% output variables
%------------------------------

if nextraout > 0
    
    extraout = cell(nt, nextraout);
    
    % Rerun at each output time
    
    for it = 1:nt
        yslice = y(it,:);
        yslice = reshape(yslice, sz);
        [extraout{it,:}] = fun(t(it), yslice, varargin{:});
    end
    
    % Permute so output arrays have time in first dimension
    
    for ivar = 1:nextraout
        ndim = ndims(extraout{1,ivar});
        varargout{ivar} = permute(cat(ndim+1, extraout{:,ivar}), [ndim+1 1:ndim]);
    end
        
end


%------------------------------
% Subfunction: Rewrite ODE 
% function so input and output 
% is vector
%------------------------------

function dydt = newfun(t, y, fun, sz, varargin)
y = reshape(y, sz);
dydt = feval(fun, t, y, varargin{:});
dydt = dydt(:);






