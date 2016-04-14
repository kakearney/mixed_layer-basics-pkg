function varargout = npz(action, varargin)
%NPZ Nutrient-phytoplankton-zooplankton biological module
%
% This module adds a simple nutrient-phytoplankton-zooplankton model to
% the mixed_layer model.  The model is based on the NPZ model described
% in Sarmiento and Gruber (2008) Chapter 4, where
%
%   dZ/dt = Z * (gammaz * g * P/(Kp + P) - lambdaz)
%   dP/dt = P * (Vmax * N/(Kn + N) - lambdap - g*Z/(Kp + P))
%   dN/dt = P * (-Vmax * N/(Kn + N) + mup * lambdap) +
%           Z * muz * ((1 - gammaz) * g * P/(Kp + P) + lambdaz)
%
% The maximum growth rate Vmax is determined by temperature and light
% limitation, following Eppley (1972) and Platt and Jasby (1976),
% respectively.  The half-saturation constants (Kn, Kp), loss rates
% (lambdap, lambdaz), and fractions of remineralization (mup, muz) are
% provided by the user.
%
% Sarmiento J, Gruber N (2006) Ocean Biogeochemical Dynamics. Princeton
% University Press 
%
% See biomodule.m for function syntax descriptions.  The following 
% fields must be present in the In structure (passed to mixed_layer as
% parameter/value pairs):
%
%   n:          n x 2 depth profile of initial nutrients, where column 1
%               gives the depth values (negative down) and column 2 holds
%               the concentrations of nutrients (mmol N m^-3)
%
%   p:          n x 2 depth profile of phytoplankton, where column 1 gives
%               the depth values (negative down) and column 2 holds the
%               concentrations of phytoplankton (mmol N m^-3)
%
%   z:          n x 2 depth profile of zooplankton, where column 1 gives
%               the depth values (negative down) and column 2 holds the
%               concentrations of zooplankton (mmol N m^-3)
%
%   kn:         half-saturation for nutrient uptake by phytoplankton (mmol
%               N m^-3)
%
%   kp:         half-saturation for phytoplankton uptake by zooplankton
%               (mmol N m^-3) 
%
%   lambdap:    loss rate for phytoplankton (s^-1)
%
%   lambdaz:    loss rate for zooplankton (s^-1)
%
%   mup:        fraction of phytoplankton loss that is remineralized
%
%   muz:        fraction of zooplankton loss that is remineralized
%
%   g:          maximum zooplankton growth rate (s^-1)
%
%   gammaz:     zooplankton assimilation efficiency (0-1)

% Copyright 2008 Kelly Kearney


nin(1) = nargin(@init) + 1;
nin(2) = nargin(@sourcesink) + 1;
nin(3) = nargin(@vertmove) + 1;

nout(1) = nargout(@init);
nout(2) = nargout(@sourcesink);
nout(3) = nargout(@vertmove);

switch action
    case 'init'
        
        out = cell(1, nout(1));
        narginchk(nin(1), nin(1));       
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
        narginchk(nin(2), nin(2));       
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
        narginchk(nin(3), nin(3));       
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};        
        
%**************************************************************************

function [bio, ismixed, bottomval, Biovars, names, diagnames] = init(In, Grd)

% Check input

p = inputParser;

p.addParamValue('n', [], @(x) size(x,2)==2);
p.addParamValue('p', [], @(x) size(x,2)==2);
p.addParamValue('z', [], @(x) size(x,2)==2);
p.addParamValue('kn', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('kp', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('lambdap', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('lambdaz', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('mup', [], @(x) isscalar(x) && isnumeric(x) && x <= 1 && x>= 0);
p.addParamValue('muz', [], @(x) isscalar(x) && isnumeric(x) && x <= 1 && x>= 0);
p.addParamValue('g', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('gammaz', [], @(x) isscalar(x) && isnumeric(x) && x <= 1 && x>= 0);

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

default = p.UsingDefaults;
if ~isempty(default)
    missing = sprintf('%s, ', default{:});
    error('Missing input parameter for npz biology: %s', missing(1:end-2));
end

In = mergestruct(p.Unmatched, p.Results);

% Two state variables, N and P, both mixed

bio(:,1) = interp1(In.n(:,1), In.n(:,2), Grd.z);
bio(:,2) = interp1(In.p(:,1), In.p(:,2), Grd.z);
bio(:,3) = interp1(In.z(:,1), In.z(:,2), Grd.z);

ismixed = true(1,3);

% Force deep-water nutrients, not phytoplankton or zooplankton

bottomval = [In.n(end,2) NaN NaN];  

% Extra variables

Biovars.kn = In.kn;
Biovars.kp = In.kp;
Biovars.mup = In.mup;
Biovars.muz = In.muz;
Biovars.lambdap = In.lambdap;
Biovars.lambdaz = In.lambdaz;
Biovars.g = In.g;
Biovars.gammaz = In.gammaz;

% Names

names = {...
	'N', 'nutrient', 'mmol N m^-3' 
	'P', 'phytoplankton', 'mmol N m^-3' 
	'Z', 'zooplankton', 'mmol N m^-3'};

% Diagnostics

diagnames = {...
	'npFlux', 'flux from N to P', 'mmol N m^-3 s^-1' 
	'pnFlux', 'flux from P to N', 'mmol N m^-3 s^-1' 
	'pzFlux', 'flux from P to Z', 'mmol N m^-3 s^-1' 
	'znFlux', 'flux from Z to N', 'mmol N m^-3 s^-1' };

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)
                                            
% Integrate using Runge-Kutta solver over the full timestep

[tout,newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss] = ...
    odewrap(@ode4, @npzode, [G.t G.t+G.dt], oldbio, [], P.T, P.par24, ...
            -G.z, G.dz, B.kn, B.kp, B.mup, B.muz, B.lambdap, B.lambdaz, ...
            B.g, B.gammaz);
                    
[newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss] = ...
    endonly(newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss);

   
% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss] = ...
    odewrap(@ode45, @npzode, [G.t G.t+G.dt], oldbio, [], P.T, P.par24, ...
            -G.z, G.dz, B.kn, B.kp, B.mup, B.muz, B.lambdap, B.lambdaz, ...
            B.g, B.gammaz);

    [newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss] = ...
        endonly(newbio,dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss);
end

% If still no good, error

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end
    
diag = [phytoUptake, phytoLoss, zpGraze, zooLoss];

%-----------------------------
% Rate of change function
%-----------------------------

function [dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss] = npzode(t, bio, temp, irr, z, dz, kn, kp, mup, muz, lambdap, lambdaz, g, gammaz)

nutrients = bio(:,1);
phyto     = bio(:,2);
zoo       = bio(:,3);
 
psRate = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz);

% Fluxes between groups

phytoUptake = psRate .* phyto;
phytoLoss = phyto .* lambdap;
zpGraze = g .* zoo .* phyto ./ (kp + phyto);
zooLoss = zoo .* lambdaz;

dbdt = zeros(size(bio));

dbdt(:,1) = mup.*phytoLoss + muz.*(zooLoss + (1-gammaz).*zpGraze) - phytoUptake; 
dbdt(:,2) = phytoUptake - phytoLoss - zpGraze;
dbdt(:,3) = gammaz.*zpGraze - zooLoss;

%**************************************************************************

function wsink = vertmove(oldbio, P, B, G)

wsink = zeros(size(oldbio));





