function newtracer = mixtracer(tracer, mixcoef, dt, dz, varargin)
%MIXTRACER Calculates diffusive mixing of a tracer in the mixed_layer model
%
% newtracer = mixtracer(tracer, mixcoef, dt, dz, sbc, bbc, bottomval,
%                       source, dissipate)
%
% This function calculates an implicit solution to integrate the diffusion
% equation over one time step.  The diffusion equation in this case is
% du/dt = d/dz(K * du/dz) - cu, where u is any tracer, K is the mixing
% coefficient for that tracer, and c is a constant dissipation term.  
%
% Input variables:
%
%   tracer:     n x 1 array of tracer concentrations at each depth at the
%               current time step (units vary based on tracer)
%
%   mixcoef:    n x 1 array, mixing coefficient for tracer at each depth
%               (m^2 s^-1) 
%
%   dt:         time increment (s)
%
%   dz:         depth increment (m)
%
% Optional input variables (passed as parameter/value pairs):
%
%   sflux:      flux of tracer across the surface interface, i.e. K * du/dz
%               at the surface, where u is tracer concentration and K is
%               the mixing coefficient (tracer unit s^-1)
%
%   bflux:      flux of tracer across the bottom interface, i.e. K * du/dz
%               along the bottom, where u is tracer concentration and K is
%               the mixing coefficient (tracer unit s^-1)
%
%   sval:       tracer value at the surface, used to force the surface grid
%               cell.
%
%   bval:       tracer value along the bottom, used to force the bottom
%               grid cell
%
%   source:     n x 1 array, source (or sink) flux of tracer at each depth
%               (tracer unit s^-1)
%
%   dissipate:  dissipation constant (s^-1)
%
% Output variables:
%
%   newtracer:  n x 1 array, tracer concentrations at next time step
%               (tracer unit)

% Copyright 2008 Kelly Kearney

%--------------------------
% Check input
%--------------------------

if ~isequal(size(tracer), size(mixcoef))
    error('Mixing coefficient array must be same length as tracer array');
end

nz = length(tracer);

% Defaults for optional parameters

A.sflux     = NaN;
A.bflux     = NaN;
A.sval      = NaN;
A.bval      = NaN;
A.source    = zeros(nz,1);
A.dissipate = zeros(nz,1);

% Parse optional parameters

pv = reshape(varargin, 2, []);

for iparam = 1:size(pv,2)
    A.(pv{1,iparam}) = pv{2,iparam};
end

%--------------------------
% Set up matrices for 
% discretized equation
%--------------------------

% Set up coefficients for tridiagonal matrix, matrix form of the equation
% a * u(j-1,t+1) + b * u(j,t+1) + c * u(j+1,t+1) = u(j,t).
% where a = -k(j)*dt/dz, b = 1 + (k(j) + k(j+1))*dt/dz + c*dt, and c =
% -k(j+1)*dt/dz

kj = mixcoef;
kjp1 = [mixcoef(2:end); NaN];

a = -kj .* (dt./dz.^2);
b = 1 + (kj + kjp1) .* (dt./dz.^2) + A.dissipate .* dt;
c = -kjp1 .* (dt./dz.^2);
d = tracer + A.source .* dt;

% Adjust coefficients for flux boundary conditions
if ~isnan(A.sflux)
    a(1) = 0;
    b(1) = 1 + kj(2) .* (dt./dz(1).^2);
    c(1) = -kj(2) .* (dt./dz(1).^2);
    d(1) = tracer(1) + A.source(1) .* dt + A.sflux .* dt;
end

if ~isnan(A.bflux)
    a(end) = -kj(end) .* (dt./dz(end).^2);
    b(end) = 1 + kj(end) .* (dt./dz(end).^2);
    c(end) = 0;
    d(end) = tracer(end) + A.source(end) .* dt - A.bflux .* dt;
end

% Adjust coefficients for value boundary conditions

if ~isnan(A.sval)
    a(1) = 0;
    b(1) = 1;
    c(1) = 0;
    d(1) = A.sval;
end

if ~isnan(A.bval)
    a(end) = 0;
    b(end) = 1;
    c(end) = 0;
    d(end) = A.bval;
end

% Create sparse tridiagonal matrix

abc = [[a(2:end); NaN] b [NaN; c(1:end-1)]];
abc = spdiags(abc, -1:1, nz, nz);

% If bottom forcing included, add additional equation,  
% u(jbottom,t+1) = bottomval, to matrix
% 
% if isnan(bottomval)
%     abc = spdiags(abc, -1:1, nz, nz);
% else
%     abc = spdiags(abc, -1:1, nz+1, nz);
%     abc(end,end) = 1;
%     d = [d; bottomval];
% end
    
%--------------------------
% Solve for u(j,t+1)
%--------------------------

newtracer = abc\d;
