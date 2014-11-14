function varargout = npzd(action, varargin)
%NPZD Nutrient-phytoplankton-zooplankton-detritus w/ upwelling
%
% [bioinit, ismixed, bottomval, Biovars, names] = biomodule('init', ...
%                                                              In, Grd);
% [newbio, diag] = biomodule('sourcesink', oldbio, meanqi, temp, z, dz, ...
%                    Biovars, t, dt);
% wsink = biomodule('vertmove', oldbio, meanqi, temp, z, dz, ...
%                    Biovars, t, dt);
%
% This module adds a nutrient-phytoplankton-zooplankton-detritus model.
% The dynamics for phytoplankton and zooplankton growth are the same as the
% NPZ model.  Egestion and mortality go to a detritus pool, which is
% remineralized at a constant rate to the nutrient pool.  Excretion goes
% directly to the nutrient pool.
%
% Input variables for 'init' mode:
%
%   In:         1 x 1 structure holding user-supplied info.  
%               
%               n:      n x 2 depth profile of initial nutrients, where
%                       column 1 gives the depth values (negative down) and
%                       column 2 holds the concentrations of nutrients
%                       (mmol N m^-3)
%
%               p:      n x 2 depth profile of phytoplankton, where column
%                       1 gives the depth values (negative down) and column
%                       2 holds the concentrations of phytoplankton (mmol N
%                       m^-3)
%
%               z:      n x 2 depth profile of zooplankton, where column
%                       1 gives the depth values (negative down) and column
%                       2 holds the concentrations of zooplankton (mmol N
%                       m^-3)
%
%               d:      n x 2 depth profile of detritus, where column
%                       1 gives the depth values (negative down) and column
%                       2 holds the concentrations of detritus (mmol N
%                       m^-3)
%
%               kn:     half-saturation for nutrient uptake by
%                       phytoplankton (mmol N m^-3) 
%
%               kp:     half-saturation for phytoplankton uptake by
%                       zooplankton (mmol N m^-3)
%
%               lambdap:loss rate for phytoplankton (s^-1)
%
%               lambdaz:loss rate for zooplankton (s^-1)
%
%               g:      maximum zooplankton growth rate (s^-1)
%
%               gammaz: zooplankton assimilation efficiency (0-1)
%
%               dsink:  detritus sinking velocity, negative (m/s)
%
%               kremin: remineralization rate constant (s^-1) 
%
%               egestz: fraction of zooplankton ingestion that is egested
%
%               upwell: upwelling speed for nutrients (m/s)
%
%   Grd:        1 x 1 structure, see mixed_layer code for fields
%
% Output variables for 'init' mode:
%
%   bioinit:    ndepth x nbsv array, initial depth profiles for all
%               biological state variables
%
%   ismixed:    1 x nbsv logical array, indicates whether each biological
%               state variable undergoes physical mixing.
%
%   bottomval:  1 x nbsv array, value used to force bottom grid cell during
%               physical mixing.  NaN indicates no forcing.
%
%   Biovars:    1 x 1 structure holding additional parameters
%
%   names:      1 x nbsv cell array of strings, names for each state
%               variable to be used as fieldname in output structure
%
%   diagnames:  1 x ndiag cell array of strings, names for each diagnostic
%               variable.  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
% Input variables for 'sourcesink' mode:
%
%   oldbio:     ndepth x nbsv array, profiles of biological state variables
%               at current time step
%
%   meanqi:     mean incoming solar radiation (W m^-2)
%
%   temp:       ndepth x 1 array, temperature profile (deg C)
%
%   z:          ndepth x 1 array, depth grid for model (m)
%
%   dz:         depth grid cell size (m)
%
%   t:          current time (s)
%
%   dt:         time step (s)
%
% Output variables for 'sourcesink' mode:
%
%   newbio:     ndepth x nbsv array, profiles of biological state variables
%               after stepping forward by one time step
%
%   diag:       ndepth x ndiag array, values of each diagnostic variable at
%               each depth
%
% Output variables for 'vertmove' mode:
%
%   wsink:      ndepth x nbsv array, velocity of vertical movement, where
%               positive values indicate movement towards the surface (m/s)

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
        error(nargchk(nargin,nin(1),nin(1)));        
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
        error(nargchk(nargin,nin(2),nin(2)));        
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
        error(nargchk(nargin,nin(3),nin(3)));        
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};
        
        
%**************************************************************************

function [bio, ismixed, bottomval, Biovars, names, diagnames] = init(In, Grd)

% Four state variables, N, P, Z, and D, all mixed

bio(:,1) = interp1(In.n(:,1), In.n(:,2), Grd.z);
bio(:,2) = interp1(In.p(:,1), In.p(:,2), Grd.z);
bio(:,3) = interp1(In.z(:,1), In.z(:,2), Grd.z);
bio(:,4) = interp1(In.d(:,1), In.d(:,2), Grd.z);

ismixed = true(1,4);

% Force deep-water nutrients, not phytoplankton, zooplankton, or detritus

% bottomval = [In.n(end,2) NaN NaN NaN];  
bottomval = nan(1,4);

% Extra variables

Biovars.kn = In.kn;
Biovars.kp = In.kp;
Biovars.lambdap = In.lambdap;
Biovars.lambdaz = In.lambdaz;
Biovars.g = In.g;
Biovars.gammaz = In.gammaz;
Biovars.dsink = In.dsink;
Biovars.kremin = In.kremin;
Biovars.egestz = In.egestz;
Biovars.upwell = In.upwell;

% Names

names = {...
	'N', 'nutrient', 		'mmol N m^-3' 
	'P', 'phytoplankton', 	'mmol N m^-3' 
	'Z', 'zooplankton', 	'mmol N m^-3',  ...
	'D', 'detritus', 		'mmol N m^-3'};
	
% Diagnostics

diagnames = {...
	'dndt', 		'dN/dt', 									'mmol N m^-3 s^-1'
 	'dpdt', 		'dP/dt', 									'mmol N m^-3 s^-1'
 	'dzdt', 		'dZ/dt', 									'mmol N m^-3 s^-1'
	'dddt', 		'dD/dt', 									'mmol N m^-3 s^-1'
	'phytoUptake', 	'uptake of N by P', 						'mmol N m^-3 s^-1'
	'phytoLoss', 	'flux of P to D due to mortality', 			'mmol N m^-3 s^-1'
	'zpGraze', 		'flux of P to Z due to grazing', 			'mmol N m^-3 s^-1'
	'zooLoss', 		'flux of Z to D due to mortality', 			'mmol N m^-3 s^-1'
	'remin', 		'flux of D to N due to remineralization', 	'mmol N m^-3 s^-1'
	'zooGrow', 		'flux of P that stays with Z as growth', 	'mmol N m^-3 s^-1'
	'zooEgest', 	'flux from P to Z to D due to egestion', 	'mmol N m^-3 s^-1'
	'zooExcrete', 	'flux from P to Z to N due to excretion', 	'mmol N m^-3 s^-1'};

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Integrate using Runge-Kutta solver over the full timestep

diag = cell(1,9);
params = {temperature, meanqi, -z, dz, Biovars.kn, Biovars.kp, ...
          Biovars.lambdap, Biovars.lambdaz, Biovars.g, Biovars.gammaz, ...
          Biovars.kremin, Biovars.egestz};

[tout,newbio, diag{:}] = odewrap(@ode4, @biochange, [t t+dt], oldbio, [], params{:});           
[newbio, diag{:}] = endonly(newbio, diag{:});
   
% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio, diag{:}] = odewrap(@ode45, @biochange, [t t+dt], oldbio, [], params{:});           
    [newbio, diag{:}] = endonly(newbio, diag{:});
end

% If still no good, error

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end

diag = cell2mat(diag);

%-----------------------------
% Rate of change function
%-----------------------------

function [dbdt, phytoUptake, phytoLoss, zpGraze, zooLoss, remin, ...
          zooGrow, zooEgest, zooExcrete] = biochange(t, bio, temp, irr, ...
          z, dz, kn, kp, lambdap, lambdaz, g, gammaz, kremin, egestz)

nutrients = bio(:,1);
phyto     = bio(:,2);
zoo       = bio(:,3);
detritus  = bio(:,4);
                                                 
% Fluxes between groups

psRate = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz);
phytoUptake = psRate .* phyto;
phytoLoss = phyto .* lambdap;
zpGraze = g .* zoo .* phyto ./ (kp + phyto);
zooLoss = zoo .* lambdaz;
remin = kremin .* detritus;
zooGrow = gammaz .* zpGraze;
zooEgest = egestz .* zpGraze;
zooExcrete = (1-gammaz-egestz) .* zpGraze;

dbdt = zeros(size(bio));

dbdt(:,1) = remin + zooExcrete - phytoUptake; 
dbdt(:,2) = phytoUptake - phytoLoss - zpGraze;
dbdt(:,3) = zooGrow - zooLoss;
dbdt(:,4) = phytoLoss + zooLoss + zooEgest - remin;

%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% Detritus sinks at a constant rate, everything else stays

wsink = zeros(size(oldbio));
wsink(:,4) = Biovars.dsink;
wsink(:,1) = Biovars.upwell;




