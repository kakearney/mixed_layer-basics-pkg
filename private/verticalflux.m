function newtracer = verticalflux(tracer, wsink, dt, dz, openbot)
%VERTICALFLUX Calculates vertical movement within mixed layer model
%
% newtracer = verticalflux(tracer, wsink, dt, dz, openbot)
%
% Calculated changes in tracer concentration due to non-mixing processes.
%
% Input variables:
%
%   tracer:     nz x 1 array, concentration of tracer
%
%   wsink:      nz x 1 array, vertical velocities (m/s)
%
%   dt:         scalar, model time step (s)
%
%   dz:         scalar, depth interval (m)
%
%   openbot:    logical scalar, true if open bottom (i.e. things sink out)
%
% Output variables:
%
%   newtracer:  nz x 1 array, new value of tracer concentrations after one
%               time step 

% Copyright 2009 Kelly Kearney

% Calculate amount of tracer leaving each depth layer

if any(tracer < 0)
    error('Negative tracer');
end

fluxout = wsink .* (dt./dz) .* tracer;

% No flux through ocean surface or floor (unless open bottom)

if fluxout(1) > 0
    fluxout(1) = 0;
end

if ~openbot && fluxout(end) < 0
    fluxout(end) = 0;
end

% Characterize flux out of each layer as going up or down

fluxdown = zeros(size(tracer));
fluxup   = zeros(size(tracer));

isup = fluxout > 0;
fluxup(isup) = fluxout(isup);
fluxdown(~isup) = -fluxout(~isup);

% Calculate total change due to loss from a layer and gain from adjacent
% layers

tracerchange = zeros(size(tracer));
tracerchange(2:end-1) = fluxup(3:end) + fluxdown(1:end-2) - abs(fluxout(2:end-1));
tracerchange(1)       = fluxup(2)                         - abs(fluxout(1));
tracerchange(end)     =                 fluxdown(end-1)   - abs(fluxout(end));

newtracer = tracer + tracerchange;
