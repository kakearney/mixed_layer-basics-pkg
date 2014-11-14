function [boygr, shear, len, gh, q2, q2l, lc_q2, lc_q2l] = mixke(ustau, vstau, ubtau, vbtau, sig, z, b1, u, v, dz, dt, gh, len, Km, lc_q2, q2, Kh, Kq, zp, zbot, q2l)
%MIXKE Calculate and mix terms related to turbulence and energy
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
%   b1:     length scale term (?)
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
%   lc_q2:  dissipation constant for turbulent kinetic energy (no units)
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
%   q2l:    (nz+1) x 1 array, turbulent kinetic energy * length scale term
%           (m^3 s^-2) 
%
% Output variables:
%
%   boygr:  (nz+1) x 1 array, buoyancy generation term
%
%   shear:  (nz+1) x 1 array, shear term
%
%   len:    (nz+1) x 1 array, turbulence length scale (m)
%
%   gh:     (nz+1) x 1 array, Richardson number (no units)
%
%   q2:     (nz+1) x 1 array, twice the turbulent kinetic energy (m^2 s^-2)
%
%   q2l:    (nz+1) x 1 array, q2 * len (m^3 s^-2)
%
%   lc_q2:  dissipation constant for turbulent kinetic energy (no units)
%
%   lc_q2l: dissipation constant for kinetic-energy-length-scale (no units)

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

cbcnst = 100;       % M&B, 6B
surfl = 2e5;        % M&B, 6B, M&B call this constant Beta

grav = 9.8;         % gravity
kappa = 0.41;       % von Karman's constant
small = 1e-9;       % pom2k small value

%------------------------------
% Surface value due to wave 
% breaking
%------------------------------

% Boundary conditions for twice the TKE ~ shear velocity squared.
% no wave breaking
% Mmntm.sbc = Mmntm.b1^(2/3)*sqrt( (Wnd.ustau(it)/Ts.Sig(1))^2 + (Wnd.vstau(it)/Ts.Sig(1))^2 );
% wave breaking, see Eq. 10 of M&B (its in the addendum)

sbc = (15.8*cbcnst)^(2/3) * sqrt((ustau/sig(1))^2 + (vstau/sig(1))^2);
  
% If wave breaking, also calculate the surface length scale as part of the
% boundary forcing, see equation 6a of M&B

l_0 = surfl*sqrt((ustau/sig(1))^2 + (vstau/sig(1))^2)/grav;

% Adjust length scale near surface according to surface waves, see eq. 5a
% of M&B.

len(1) = kappa * l_0;
isnearsurf = z < 0 & z > -0.5;
len(isnearsurf) = max(len(isnearsurf), kappa * l_0);

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

if isscalar(dz)
    q2 = mixtracer(q2, [0; Kqt], dt, dz, 'sval', sbc, 'bval', bbc, ...
                  'source', 2*ss*(dz./dt), 'dissipate', 2*lc_q2); 
else
    dzextended = [dz; dz(end)];
    q2 = mixtracer(q2, [0; Kqt], dt, dzextended, 'sval', sbc, 'bval', bbc, ...
                  'source', 2.*ss.*(dzextended./dt), 'dissipate', 2*lc_q2); 
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

wallprox(2:nz) = 1 + e2*(len(2:nz).*(1./(abs(zp(2:nz))) + ...
                         1./abs(zp(2:nz)-zbot)) ./ (abs(zbot)*kappa)).^2;
wallprox(1)   = 0;
wallprox(end) = 0;

lc_q2l = lc_q2 .* wallprox;

if isscalar(dz)
    q2l = mixtracer(q2l, [0; Kqt], dt, dz, 'sval', sbc, 'bval', bbc, ...
                'source', ss*(dz./dt), 'dissipate', lc_q2l);
else
    q2l = mixtracer(q2l, [0; Kqt], dt, dzextended, 'sval', sbc, 'bval', bbc, ...
                'source', ss.*(dzextended./dt), 'dissipate', lc_q2l);
end

%------------------------------
% Adjust for small values
%------------------------------

issmall = q2 < small;
issmall(1) = false;
issmall(nz+1) = false;
q2(issmall) = small;

issmall = q2l < small;
issmall(1) = false;
issmall(nz+1) = false;
q2l(issmall) = small;

%------------------------------
% Calculate length scales and 
% Reynold's numbers based
% on new q2 and q2l
%------------------------------

% Derive for all but top and bottom (length scale = 0 for these)

len(2:nz) = q2l(2:nz)./q2(2:nz);
gh(2:nz) = (len(2:nz).^2./q2(2:nz)).*boygr(2:nz);

% Length scale limitation

islim = (len.^2 .* boygr) < (-0.281 .* q2);
islim(1) = false;
islim(nz+1) = false;
q2l(islim) = q2(islim) .* sqrt(-0.281*q2(islim)./(boygr(islim) - small));
len(islim) = q2l(islim)./q2(islim);
gh(islim) = (len(islim).^2./q2(islim)).*boygr(islim);

% not too sure how this would occur - I guess static instabilities can be
% handles automoatically - response is elevated mixing.  Will need to try
% code with and without static instability adjustment.

gh(gh > 0.028) = 0.028;
