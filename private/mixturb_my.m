function [boygr, shear, len, gh, q2, q2l, lc_q2, lc_q2l] = ...
         mixturb(ustau, vstau, ubtau, vbtau, sig, z, b1, u, v, dz, dt, ...
               gh, len, Km, lc_q2, q2, Kh, Kq, zp, zbot, q2l)
%MIXTURB_MY Calculate and mix terms related to turbulence and energy
%
% [boygr, shear, len, gh, q2, q2l, lc_q2, lc_q2l] = mixke(ustau, vstau, ...
% ubtau, vbtau, sig, z, b1, u, v, dz, dt, gh, len, Km, lc_q2, q2, Kh, ...
% Kq, zp, zbot, q2l)
%
% This function calculates several terms related to the turbulence closure
% scheme used in this model.  It calculates bouyancy and shear values based
% on the wind stress and water density at the current time step, then uses
% this to mix the turbulent kinetic energy and kinetic-energy-length-scale
% terms (q2 and q2l) and to calculate new length scales and Richardson
% numbers for each depth.
%
% Input variables:
%
%   ustau:  surface east-west wind stress (N m^-2)
%
%   vstau:  surface north-south wind stress (N m^-2)
%
%   ubtau:  bottom east-west stress (N m^-2)
%
%   vbtau:  bottom north-south stress (N m^-2)
%
%   sig:    nz x 1 array, seawater density (kg m^-3)
%
%   z:      nz x 1 array, depth of each cell (m)
%
%   b1:     constant (unitless)
%
%   u:      nz x 1 array, east-west water velocity (m s^-1)
%
%   v:      nz x 1 array, north-south water velocity (m s^-1)
%
%   dz:     depth increment (m)
%
%   dt:     time increment (s)
%
%   gh:     (nz+1) x 1 array, Richardson number (no units)
%
%   len:    (nz+1) x 1 array, turbulence length scale (m)
%
%   Km:     (nz+1) x 1 array, vertical kinematic viscosity, i.e. momentum
%           mixing coefficient (m^2 s^-1) 
%
%   lc_q2:  dissipation constant for turbulent kinetic energy (s^-1)
%
%   q2:     (nz+1) x 1 array, twice the turbulent kinetic energy (m^2 s^-2)
%
%   Kh:     tracer mixing coefficient (m^2 s^-1)
%
%   Kq:     turbulence mixing coefficient (m^2 s^-1)
%
%   zp:     (nz+1) x 1 array, depth of edges of grid cells (m)
%
%   zbot:   maximum depth (m)
%
%   q2l:    (nz+1) x 1 array, twice the turbulent kinetic energy * length
%           scale (m^3 s^-2)
%
% Output variables:
%
%   boygr:  (nz+1) x 1 array, buoyancy generation term (s^-2, multiplied by
%           Kh to give a rate of change of turbulent energy in m^2 s^-3)
%
%   shear:  (nz+1) x 1 array, rate of change in turbulent energy due to
%           velocity shear (m^2 s^-3)
%
%   len:    (nz+1) x 1 array, turbulence length scale (m)
%
%   gh:     (nz+1) x 1 array, Richardson number (no units)
%
%   q2:     (nz+1) x 1 array, twice the turbulent kinetic energy (m^2 s^-2)
%
%   q2l:    (nz+1) x 1 array, q2 * len (m^3 s^-2)
%
%   lc_q2:  dissipation constant for turbulent kinetic energy (s^-1)
%
%   lc_q2l: dissipation constant for kinetic-energy-length-scale (s^-1)

% Copyright 2008 Kelly Kearney

%------------------------------
% Setup parameters and 
% constants
%------------------------------
    
nz = size(sig,1);

% Zero out some arrays

ss       = zeros(nz+1, 1);
wallprox = zeros(nz+1, 1);
Kqt      = zeros(nz,   1);

% Constants used in Mmntm.q2l equation

e1 = 1.8;
e2 = 1.33;
e3 = 1.0;

% constants for wave breaking, see Mellor and Blumberg, JPO, March 2004,
% pp. 693-698.

alpha_cb = 100;       % M&B, 6B, recommended value for constant relating
                      % surface q2 to the shear velocity.
beta = 2e5;           % M&B, 6B, recommended value for constant used to
                      % determine the turbulent length-scale near the ocean
                      % surface.

grav = 9.8;         % gravity
kappa = 0.41;       % von Karman's constant
small = 1e-6;       % pom2k small value

%------------------------------
% Surface value due to wave 
% breaking
%------------------------------

% Surface boundary conditions q2 via Eq. 10 of Mellor and Blumberg (2004)
sbc = (15.8*alpha_cb)^(2/3) * sqrt((ustau/sig(1))^2 + (vstau/sig(1))^2);
sbc = max(sbc,1e-4);
  
% Surface length scale boundary condition from equation 6a of Mellor and
% Blumberg (2004)
z_w = beta*sqrt((ustau/sig(1))^2 + (vstau/sig(1))^2)/grav;
% CAS: set minimum surface length scale
z_w = max(z_w,0.02);
len(1) = kappa * z_w;
gh(1) = 0;

%------------------------------
% Bottom value (due to bottom
% stress)
%------------------------------

bbc = b1^(2/3) * sqrt((ubtau/sig(end))^2 + (vbtau/sig(end))^2);
gh(end-1) = 0;

%------------------------------
% Source/sink due to shear and
% buoyancy
%------------------------------

boygr = zeros(size(gh));
shear = zeros(size(gh));

% Calculate the changes in TKE due to shear and buoyancy

if isscalar(dz)
    boygr(2:end-1) = (grav./mean(sig)).*(sig(1:end-1) - sig(2:end))./dz;
    shear(2:end-1) = Km(2:end-1) .* (((u(1:end-1) - u(2:end))./dz).^2 + ...
                                     ((v(1:end-1) - v(2:end))./dz).^2);
else
    boygr(2:end-1) = (grav./mean(sig)).*(sig(1:end-1) - sig(2:end))./dz(2:end);
    shear(2:end-1) = Km(2:end-1) .* (((u(1:end-1) - u(2:end))./dz(2:end)).^2 + ...
                                     ((v(1:end-1) - v(2:end))./dz(2:end)).^2);
end

% Net production, shear term is always positive, boygr becomes more
% negative as stratification strengthens, factor of 2 because MY scheme
% tracks 2 times the TKE.

ss(2:nz) = shear(2:nz) + Kh(2:nz) .* boygr(2:nz);

%------------------------------
% Dissipation term
%------------------------------

% Coefficient for the linearized dissipation term.  Dissipation is
% proportional to q^3/l, the expression below is q/l for the present time
% step.  Mmntm.q2 is then solved for in the diffusion equation.  This is
% possible because dissipation is always negative, meaning the tri-diagonal
% matrix in the solution will stay diagonally dominant.  
%
% EZER (2000), MELLOR (2001) suggest a few possible corrections for this

lc_q2(1:nz) = sqrt(q2(1:nz)) ./ (b1.*len(1:nz) + small);

%------------------------------
% Mixing coefficient
%------------------------------

% Average Kq's to define a Kq variable at the grid center for the
% Thomas Algorithm. 

Kqt(1:nz) = (Kq(1:nz) + Kq(2:(nz+1))) ./ 2;

%------------------------------
% Mix q^2
%------------------------------

% CAS: fixed units issue with source/sink term, note that the factor of 2 
% for the source/sink and dissipation terms is because q^2 is 2*the 
% turbulent kinetic energy
if isscalar(dz)
    q2 = mixtracer(q2, [0; Kqt], dt, dz, 'sval', sbc, 'bval', bbc, ...
                   'source',2*ss, 'dissipate', 2*lc_q2);
else
    dzextended = [dz; dz(end)];
    q2 = mixtracer(q2, [0; Kqt], dt, dzextended, 'sval', sbc, 'bval', bbc, ...
                  'source', 2.*ss, 'dissipate', 2*lc_q2); 
end

%------------------------------
% Mix q^2*l
%------------------------------

% Boundary values
sbc = sbc * len(1);
bbc = 0;

% Source/sink from bouyancy and shear

ss(2:nz) = e1*len(2:nz).*(shear(2:nz) + e3*Kh(2:nz).*boygr(2:nz));

% If it is a large length scale and you are close to the boundary, then
% augment the dissipation term for q2l so that l will decline.

% CAS: adjustment to wall proximity function for consistency with
% Williams (2006).  
wallprox(2:nz) = 1 + e2*(len(2:nz)/kappa .* ...
                    (1./abs(zp(2:nz)) + 1./abs(zp(2:nz)-zbot)) ).^2;
wallprox(1)   = 0;
wallprox(end) = 0;

% The disipaption of q2l is calculated as wallprox*q/(b1*len) * (q2*len),
% where the first term is calculated from conditions at the last time step.
% Note that lc_q2 - q/(b1*len)
lc_q2l = lc_q2 .* wallprox;

% CAS: fixed units issue with source/sink term
if isscalar(dz)
    q2l = mixtracer(q2l, [0; Kqt], dt, dz, 'sval', sbc, 'bval', bbc, ...
                  'source',ss,'dissipate',lc_q2l);
else
    q2l = mixtracer(q2l, [0; Kqt], dt, dzextended, 'sval', sbc, 'bval', bbc, ...
                'source', ss, 'dissipate', lc_q2l);
end

%-------------------------------------------------------------------------
% Adjust for small values in the interior, do not change boundary values
%-------------------------------------------------------------------------

aa = (q2 < small) | (q2l < small);
q2(aa) = small;
q2l(aa) = 0.1*150*small;

%------------------------------
% Calculate length scales and 
% Reynold's numbers based
% on new q2 and q2l
%------------------------------

% Derive len for all but top and bottom (length scale = 0 for these)

len(2:nz) = q2l(2:nz)./q2(2:nz);

% CAS: calculate Richardson term - note that water has been treated as
% incompressible for this calculation.
gh(2:nz) = (len(2:nz).^2./q2(2:nz)).*boygr(2:nz);

% Length scale limitation in highly stratified flows see Galperin et al. 
% (1988), eq. (22)
islim = gh < -0.28;
islim(1) = false;
islim(nz+1) = false;
% Reduce the length scale such that gh = 0.28, noting that:
% q2l = q2*sqrt(gh*q2/boygr)
q2l(islim) = q2(islim) .* sqrt(-0.28*q2(islim)./(boygr(islim) - small));
len(islim) = q2l(islim)./q2(islim);

% Calculate gh

gh(islim) = (len(islim).^2./q2(islim)).*boygr(islim);

% avoids SH singularity in transient highly unstable situations
gh(gh > 0.028) = 0.028;

% Make the length scale near the surface consistent with the wave-driven
% boundary condition.  This ensures a consistent length scale near the
% surface in cases where strong dissipation of l by the wall proximity
% function over-rides the BC.
isnearsurf = z < 0 & z > -kappa*z_w;
len(isnearsurf) = max(len(isnearsurf), kappa * z_w);
