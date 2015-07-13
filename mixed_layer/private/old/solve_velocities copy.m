function [newU newV] = solve_velocities(oldvels, mixcoef, dt, dz, cor, ...
                                        pgx, pgy, varargin)
% solve_velocities steps the equations of motion forward to solve for new
% velocity values
%
% newvel = solve_velocities(oldvel, mixcoef, dt, dz,cor,pgx,pgy, ...
%                           sbc,bbc,bottomval,source, dissipate)
%
% The routine uses a semi-implicit Crank-Nicolson scheme to solve:
%
% du/dt = fv + d/dz(K*du/dz) - pgx - eps*u
% dv/dt = -fu + d/dz(K*dv/dz) - pgy - eps*v
%
% Where u and v are the velocities, f is the coriolis parameter, K is a
% mixing coefficient, pgx and pgy are pressure gradients in the x and y
% directions respectively, and eps is a dissipation term that mimics the
% horizontal divergence.  
%
% Input variables (n = number of vertical layers)
%
%   oldvel:     2n x 1 array of velocities at each depth at the current 
%               time step, first n entries are u values, next n entries are
%               v values (m s^-1)
%
%   mixcoef:    n x 1 array, mixing coefficient for tracer at each depth
%               (m^2 s^-1)
%
%   dt:         time increment (s)
%
%   dz:         depth increment (m)
%
%   cor:        coriolis parameter (s^-1)
%
%   pgx:        pressure acceleration in the x direction (m s^-2)
%
%   pgy:        pressure acceleration in the y direction (m s^-2)
%
% Optional input variables (passed as parameter/value pairs):
%
%   sflux:      "velocity flux" i.e., momentum flux/density, into the
%               top depth bin = 1/dz(1) * K * du/dz (m s^-2)
%
%   bflux:      "velocity flux" i.e., momentum flux/density, into the top
%               depth bin = 1/dz(1) * K * du/dz (m s^-2)
%
%   sval:       velocity value at the surface, used to force the surface
%               grid cell (m/s)
%
%   bval:       velocity value along the bottom, used to force the bottom
%               grid cell (m/s)
%
%   source:     n x 1 array, source (or sink) flux of tracer at each grid 
%               cell (m s^-2)
%
%   dissipate:  dissipation constant (s^-1)
%
% Output variables:
%
%   newvel:     n x 1 array, tracer concentrations at next time step
%               (tracer unit)
%
% Local variables:
%
%   theta:      Numerical parameter determining the relative weighting of
%               future versus present conditions in the numerical
%               calculation.  A value of 0 is an explicit scheme, a value
%               of 1 is a fully implicit scheme, the default value of 0.5
%               is the Crank-Nicolson scheme, which offers a robust blend
%               of stability and accuracy.
%
%   dtheta:     Added this parameter to control the relative weighting of
%               future versus present conditions for the diffusive step.
%               
%

% Copyright 2008 Kelly Kearney

%--------------------------
% Check input
%--------------------------

if ~isequal(size(oldvels,1), 2*size(mixcoef,1))
    error('Velocity array must be twice the length as tracer array');
end

%------------------------------
% Define some local variables
%------------------------------

theta = 0.5;
thetad = 1.0;

nz = length(oldvels)/2;

% Defaults for optional parameters

A.sflux_u     = NaN;
A.sflux_v     = NaN;
A.bflux_u     = NaN;
A.bflux_v     = NaN;
A.sval_u      = NaN;
A.sval_v      = NaN;
A.bval_u      = NaN;
A.bval_v      = NaN;
A.source      = zeros(2*nz,1);
A.dissipate = 0;

% Parse optional parameters

pv = reshape(varargin, 2, []);

for iparam = 1:size(pv,2)
    A.(pv{1,iparam}) = pv{2,iparam};
end

%--------------------------
% Set up matrices for 
% discretized equation
%--------------------------

%  Set up the coefficients for the sparse matrix used to solve the momentum
%  equation.  The matrix consists of a tri-diagonal portion containing
%  diffusive and dissipative terms for first u and then v.  The coriolis 
%  accelerations produce terms displaced from the main diagonal by nz 
%  spaces due to the u,v cross-terms they contain. The matrix form of the 
%  equation for the first nz rows is:
%  a * u(j-1,t+1) + b * u(j,t+1) + c * u(j+1,t+1) + d1 * v(j+1,t+1) = e
%
%  For the next nz rows we have:
%  a * v(j-1,t+1) + b * v(j,t+1) + c * v(j+1,t+1) - d2 * u(j+1,t+1) = e
%
%  The coefficients a, b, c, and d are defined below and they are arranged
%  to form the sparse matric abcd.  e is a 2n x 1 column vector whose
%  entries depend only on known values at the present time step.  New
%  velocity values are found via:
%
%  newvels = abcd\e
%

kj = [mixcoef; mixcoef];
kjp1 = [mixcoef(2:end); NaN; mixcoef(2:end); NaN];

a = -thetad * kj .* (dt/dz^2);
b = 1 + thetad * (kj + kjp1) .* (dt/dz^2) + theta * A.dissipate * dt;
c = -thetad * kjp1 * (dt/dz^2);
d_u = -theta * cor * dt * ones(size(a));
d_l = theta * cor * dt * ones(size(a));

% construct rhs from old velocity information
velj = oldvels;
veljm1 = [NaN; oldvels(1:nz-1); NaN; oldvels(nz+1:end-1)];
veljp1 = [oldvels(2:nz); NaN; oldvels(nz+2:end); NaN];
e = velj + (1-thetad)*(kj.*veljm1 - (kj + kjp1).*velj + kjp1.*veljp1) .* ...
    (dt/dz^2) - (1-theta) * A.dissipate .* velj * dt + A.source * dt;
e(1:nz) = e(1:nz) + (1-theta) * cor * velj((nz+1):2*nz) * dt - pgx * dt;
e(nz+1:2*nz) = e(nz+1:2*nz) - (1-theta) * cor * velj(1:nz) * dt - pgy * dt;

% Adjust coefficients for flux boundary conditions, 1 is the top boundary
% for u, nz+1 is the top boundary for v.

if ~isnan(A.sflux_u)
    a(1) = 0;
    b(1) = 1 + thetad * kj(2) * (dt/dz^2) + theta * A.dissipate * dt;
    c(1) = -thetad * kj(2) * (dt/dz^2);
    e(1) = velj(1) + (1-thetad) * (-kj(2)*velj(1) + kj(2)*velj(2)) * ... 
           (dt/dz^2) - (1-theta) * A.dissipate * velj(1) * dt + ...
           A.source(1) * dt + (1-theta) * cor * velj(nz+1) * dt - ...
           pgx * dt + A.sflux_u * dt;
end

if ~isnan(A.sflux_v)
    a(nz+1) = 0;
    b(nz+1) = 1 + thetad * kj(2) * (dt/dz^2) + theta * A.dissipate * dt;
    c(nz+1) = -thetad * kj(2) * (dt/dz^2);
    e(nz+1) = velj(nz+1) + (1-thetad) * (-kj(nz+2)*velj(nz+1) + ...
              kj(nz+2)*velj(nz+2)) * (dt/dz^2) - (1-theta) * ...
              A.dissipate * velj(nz+1) * dt + A.source(nz+1) * dt - ...
              (1-theta) * cor * velj(1) * dt - pgy * dt + A.sflux_v * dt;
end

if ~isnan(A.bflux_u)
    a(nz) = -thetad * kj(nz) * (dt/dz^2);
    b(nz) = 1 + thetad * kj(nz) * (dt/dz^2) + theta * A.dissipate * dt;
    c(nz) = 0;
    e(nz) = velj(nz) + (1-thetad) * (-kj(nz)*velj(nz) + ...
            kj(nz)*velj(nz-1)) * (dt/dz^2) - (1-theta) * ...
            A.dissipate * velj(nz) * dt + A.source(nz) * dt + ...
            (1-theta) * cor * velj(end) * dt - pgx * dt - A.bflux_u * dt;
end

if ~isnan(A.bflux_v)       
    a(end) = -thetad * kj(nz) * (dt/dz^2);
    b(end) = 1 + thetad * kj(nz) * (dt/dz^2) + theta * A.dissipate * dt;
    c(end) = 0;
    e(end) = velj(end) + (1-thetad) * (-kj(nz)*velj(end) + ...
             kj(nz)*velj(end-1)) * (dt/dz^2) - (1-theta) * ...
             A.dissipate * velj(end) * dt + A.source(end) * dt - ...
             (1-theta) * cor * velj(nz) * dt - pgy * dt - A.bflux_v * dt;
end

% Adjust coefficients for value boundary conditions

if ~isnan(A.sval_u)
    a(1) = 0;
    b(1) = 1;
    c(1) = 0;
    d_u(nz+1) = 0;  % spdiags takes from the lower part for super-diags
    d_l(1) = 0;     
    e(1) = A.sval_u;
end

if ~isnan(A.sval_v)
    a(nz+1) = 0;
    b(nz+1) = 1;
    c(nz+1) = 0;
    d_l(nz+1) = 0;  % spdiags takes from the upper part for sub-diags
    e(nz+1) = A.sval_v;
end

if ~isnan(A.bval_u)
    a(nz) = 0;
    b(nz) = 1;
    c(nz) = 0;
    d_u(end) = 0;               
    e(nz) = A.bval_u;
end

if ~isnan(A.bval_v)
    a(end) = 0;
    b(end) = 1;
    c(end) = 0;
    d_l(nz) = 0;
    e(end) = A.bval_v;
end

% Create sparse tridiagonal matrix

abcd = [d_l [a(2:end); NaN] b [NaN; c(1:end-1)] d_u];
abcd = spdiags(abcd, [-nz -1 0 1 nz], 2*nz, 2*nz);

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

newvels = abcd\e;
newU = newvels(1:nz);
newV = newvels(nz+1:end);
