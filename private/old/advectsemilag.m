function psinext = advectsemilag(x, psi, u, t, dt, bound)
%ADVECTSEMILAG Semi-Lagrangian advection
%
% psinew = advectsemilag(x, psi, u, t, dt)
%
% This function calculates tracer concentration following a single time
% step with semi-Lagrnagian advection in one dimension, following the
% method in Numerical Methods for Wave Equations in Geophysical Fluid
% Dynamics by Dale R. Durran.
%
% Input variables:
%
%   x:      n x 1 array, x coordinate of grid cells
%
%   psi:    n x 1, tracer value at current time in each grid cell
%
%   u:      Advection speed.  This can be either a scalar, indicating a
%           constant velocity over space and time, or a function handle to
%           a function f(t,x) that accepts and returns scalars.
%
%   t:      1 x 1, current time
%
%   dt:     1 x 1, time step
%
%   bound:  how to deal with boundaries
%           'zero':     set value for departure points outside of grid to 0
%           'nearest':  set value of departure points outside of grid to
%                       value of nearest grid cell
%           'periodic': assume spatial domain is periodic, i.e. tracer
%                       flowing out one end reenters at other
%
% Output variables:
%
%   psinew: 1 x n, tracer value at next time step    

% Copyright 2009 Kelly Kearney

%-------------------------
% Check input
%-------------------------

if isa(u, 'function_handle')
    constant = false;
else
    constant = true;
end

if nargin < 6
    bound = 'zero';
end

%-------------------------
% Calculate departure 
% point for each grid 
% point
%-------------------------

if constant
    xdepart = x - u.*dt;
else
    nx = length(x);
    u1 = u(t, x);
    xstar = x - u1.*dt/2;
    u2 = u(t+dt/2, xstar);
    xdepart = x - u2.*dt;
end

% Interpolate psi value at each departure point
% TODO: right now assuming no advection across boundaries by setting value
% of departure points outside grid to 0... does this make sense?

psinext = interp1(x, psi, xdepart, 'pchip', NaN);

% If departure point is outside grid, use nearest psi value

isout = isnan(psinext);

if any(isout)
    switch bound
        case 'zero'
            psinext(isout) = 0;
        case 'nearest'
            
            fromleft  = isout & (xdepart < min(x) | (isnan(xdepart) & u1 > 0));
            fromright = isout & (xdepart > max(x) | (isnan(xdepart) & u1 < 0));
            [ileft, ileft] = min(x);
            [iright, iright] = max(x);
            psinext(fromleft) = psi(ileft);
            psinext(fromright) = psi(iright);
            
%             if isout(end) & u1(end) > 0 % upwell through bottom
%                 psinext(end) = psi(end);
%             elseif isout(1) & u1(1) < 0 % downwell through surface
%                 psinext(1) = psi(1);
%             end
%             psiextrap = interp1(x, psi, xdepart, 'nearest', 'extrap');
%             psinext(isout) = psiextrap(isout);
        case 'periodic'
            xrange = diff(minmax(x));
            xdepart2 = mod(xdepart - min(x), xrange) + min(x);
            psinext = interp1(x, psi, xdepart2, 'pchip');
        case 'periodic2'
            xrange = diff(minmax(x));
            xdepart2 = mod(xdepart - min(x), xrange) + min(x);
            xdepart2(xdepart2 == min(x)) = max(x);
            psinext = interp1(x, psi, xdepart2, 'pchip');
        otherwise
            error('Unrecognized boundary method');
    end
end

