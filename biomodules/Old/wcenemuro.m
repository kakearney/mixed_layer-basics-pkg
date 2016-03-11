function varargout = wcenemuro(action, varargin)
%WCENEMURO Water column ecosystem with NEMURO biological module
%
% [bioinit, ismixed, bottomval, Biovars, names] = ...
%    wcenemuro('init', In, Grd);
% [newbio, diag] = ...
%    wcenemuro('sourcesink', oldbio, meanqi, temp, z, dz, Biovars, t, dt);
% wsink = ...
%    wcenemuro('vertmove', oldbio, meanqi, temp, z, dz, Biovars, t, dt);
%
% This module runs a full end-to-end ecosystem model.  The lower trophic
% levels (nutrient cycling, primary production, and lowest-level grazers)
% are based on the NEMURO model, while upper trophic level interactions
% are based on the Ecopath model,  The biological state variables include
% both planktonic and nektonic functional groups; planktonic groups are
% confined to a specific layer of the water column model, while nektonic
% groups can freely move (and feed) throughout the water column (for
% computational purposes this also includes non-water-dwelling predators
% like birds).
%
% Input variables for 'init' mode:
%
%   In:         1 x 1 structure holding user-supplied info.  
%
%               Ewein :     Ewe input structure (see
%                           ecopathlite.m).References to ngroup below
%                           refers to the number of functional groups in
%                           this model.
%
%               nekidx:     vector of indices corresponding to the groups
%                           in the Ewe model that will be nektonic in the
%                           water column model (i.e. can feed on the
%                           entire water column).
%
%               weightfun:  function handle to function of the format w =
%                           f(z) describing the distribution of profile
%                           shape of plankton concentration at time 0. 
%                           The values will be normalized to distribute
%                           the steady-state biomass throughout the water
%                           column appropriately. z represents depth in
%                           meters (positive down)
%
%               switching:  ngroup x 1, switching parameter for each
%                           predator. 0 = no switching, < 1 = switches
%                           only when prey is very rare, > 1 switches
%                           quickly (no units, 0-2)
%
%               bioinit:    nbz x 11 array, initial biomass values for
%                           all NEMURO-derived variables, i.e. PS, PL,
%                           ZS, ZL, ZP, NO3, NH4, PON, DON, SiOH4, and
%                           Opal. If NaNs are found rather than a value,
%                           then the default initial value from the
%                           original NEMURO fortran code is used
%                           instead.(mol/l)
%
%               bioinitz:   nbz x 1 array, vector of depth value
%                           corresponding to the rows of bioinit (m,
%                           positive down (yes, that's confusing, differs
%                           from the physical stuff due to the way nemuro
%                           was coded originally))
%
%               nemidx:     1 x 11 vector, indices of Ewein groups that
%                           correspond to each NEMURO-derived variable. 
%                           A NaN indicates that there is no functional
%                           group corresponding to this variable.
%
%               diagnostic: n x 1 cell array of strings, where each
%                           string is a regular expression.  Diagnostic
%                           fields matching any of these regular
%                           expressions will be included in the output
%                           file.  Diagnostic fields are of the format
%                           axxyy, where a is one of the strings dbdt,
%                           primprod, respiration, decomp, nitrify, pred,
%                           egestion, excretion, or mortality, and xx and
%                           yy are two-digit indices of the source and
%                           sink functional group, respectively.
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
%   bottomval:  1 x nbsv array, value used to force bottom layer during
%               mixing.  NaN indicates no forcing.
%
%   Biovars:    1 x 1 structure holding additional parameters
%
%   names:      1 x nbsv cell array of strings, names for each state
%               variable,  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
%   diagnames:  1 x ndiag cell array of strings, names for each diagnostic
%               variable.  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
% Input variables for 'sourcesink' and 'vertmove' mode:
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


%-------------------------
% Set up biological state
% variable initial 
% conditions and names,
% and additional 
% parameters
%-------------------------

[bio, names, Biovars] = wcenemurosetup(In, Grd);

%-------------------------
% Mixing and deep water 
% forcing
%-------------------------

% Don't mix nektonic groups

ismixed = ~Biovars.isnekton;

% No deep water forcing

bottomval = nan(Biovars.ngroup,1);

%-------------------------
% Diagnostics
%-------------------------

% Source and sink fluxes

diagtokeep = {'dbdt', 'primprod', 'respiration', 'decomp', 'nitrify', ...
              'pred', 'egestion', 'excretion', 'mortality'};
          
namesplus = [names; {'--', 'critter silica', '--'}];

count = 0;
for id = 1:length(diagtokeep)
    if strcmp(diagtokeep{id}, 'dbdt')
        for ig = 1:Biovars.ngroup
            count = count + 1;
            diagnames{count,1} = sprintf('%s%02d', diagtokeep{id}, ig);
            diagnames{count,2} = sprintf('%s, %s', diagtokeep{id}, names{ig,2});
            diagnames{count,3} = 'mmol N m^-2 s^-1';
        end
    else
        for ig2 = 1:(Biovars.ngroup+1)
            for ig1 = 1:(Biovars.ngroup+1)
                count = count + 1;
                diagnames{count,1} = sprintf('%s%02d%02d', diagtokeep{id}, ig1, ig2);
                diagnames{count,2} = sprintf('%s, %s to %s', diagtokeep{id}, namesplus{ig1,2}, namesplus{ig2,2});
                diagnames{count,3} = 'mmol N m^-2 s^-1';
            end
        end
    end
end

% Flux matrices are usually sparse, and keeping all of them leads to huge
% files, so only keep the ones the user indicates 

if ~isfield(In, 'diagnostic')
    In.diagnostic = {};
end

keep = false(size(diagnames,1),1);
for ikeep = 1:length(In.diagnostic)
    keep = keep | regexpfound(diagnames(:,1), In.diagnostic(ikeep));
end
Biovars.keepdiag = keep;
Biovars.diagnames = diagnames; % Debugging
diagnames = diagnames(keep,:);


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Set up input parameters for ode function

% WceParam = rmfield(Biovars, {'ngroup', 'wsink', 'loss', 'keepdiag', 'forceIngest'});
% if isfield(Biovars, 'time')
%     WceParam = rmfield(WceParam, 'time');
% end
WceParam.z = -z;
WceParam.dz = dz;
WceParam.irr = meanqi;
WceParam.temp = temperature;
WceParam.other = Biovars.loss;

WceParam = Biovars;
WceParam.z = -z;
WceParam.dz = dz;
WceParam.irr = meanqi;
WceParam.temp = temperature;
WceParam.other = Biovars.loss;
      
% Forcing (if used)
% 
% if ~isempty(Biovars.forceIngest)
%     it = Biovars.time == t;
%     nf = size(Biovars.forceIngest,1);
%     WceParam.forceIngest = zeros(nf, length(z)+2);
%     for ig = 1:nf
%         WceParam.forceIngest(ig,1) = Biovars.forceIngest{ig,1};
%         WceParam.forceIngest(ig,2) = Biovars.forceIngest{ig,2};
%         WceParam.forceIngest(ig,3:end) = Biovars.forceIngest{ig,3}(:,it);
%     end    
% else
%     WceParam.forceIngest = [];
% end
    
% Convert bio from concentration/m^3 to conc/m^2

oldbio = bsxfun(@times, oldbio, dz);
      
% Integrate using Runge-Kutta solver over the full timestep

odefun = @(t,b,P) wcenemuroode(t, b, P); % specify 3 rather than 4 inputs
[tout, newbio, db, WceOut] = odewrap(@ode4, odefun, [t t+dt], oldbio, [], WceParam);           

[newbio, diag{1}] = endonly(newbio, db);
diag = [diag; struct2cell(WceOut(1))]';


% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

shouldbenan = false(size(oldbio));
shouldbenan(2:end, Biovars.isnekton) = true;

Opt = odeset('nonnegative', 1:numel(oldbio));
if any((isnan(newbio(:)) & ~shouldbenan(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout, newbio, db, WceOut] = odewrap(@ode45, odefun, [t t+dt], oldbio, Opt, WceParam);  
    diag = cell(1,1);
    [newbio, diag{1}] = endonly(newbio, db);
    if any(newbio(shouldbenan))
        error('Numbers where should be NaN placeholder');
    else
        newbio(shouldbenan) = NaN;
    end
    diag = [diag; struct2cell(WceOut(1))]';
end

% If still no good, error

if any((isnan(newbio(:)) & ~shouldbenan(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end
   
% Reformat diagnostics

for id = 1:length(diag)
    if ndims(diag{id}) == 3
        diag{id} = permute(diag{id}, [3 1 2]);
        diag{id} = reshape(diag{id}, size(diag{id},1), []);
    end
end
diag = cat(2, diag{:});
diag = diag(:, Biovars.keepdiag);

% Convert back to m^-3

newbio = bsxfun(@rdivide, newbio, dz);


%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% PON and Opal settle

ponidx = Biovars.nemidx(8);
opalidx = Biovars.nemidx(11);

biosettle = zeros(1, Biovars.ngroup);
biosettle(1,ponidx)  = -Biovars.setVPON; 
biosettle(1,opalidx) = -Biovars.setVOpal;

wsink = ones(size(z)) * biosettle;
