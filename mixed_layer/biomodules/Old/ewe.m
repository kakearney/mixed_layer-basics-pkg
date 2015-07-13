function varargout = ewe(action, varargin)
%EWE Adaptation of EwE model
%
% [bioinit, ismixed, bottomval, Biovars, names] = biomodule('init', ...
%                                                              In, Grd);
% [newbio, diag] = biomodule('sourcesink', oldbio, meanqi, temp, z, dz, ...
%                    Biovars, t, dt);
% wsink = biomodule('vertmove', oldbio, meanqi, temp, z, dz, ...
%                    Biovars, t, dt);
%
% This module adapts the Ecopath with Ecosim (EwE) model so that the
% primary production included in the food web responds to specific
% environmental factors.  It also breaks the functional groups into
% planktonic and nektonic groups, where planktonic groups are confined
% to a specific layer of the water column model, while nektonic groups can
% freely move (and feed) throughout the water column (for computational
% purposes this also includes non-water-dwelling predators like birds).
%
% Input variables for 'init' mode:
%
%   In:         1 x 1 structure holding user-supplied info.  
%               
%               !!! LIST USER INPUT VARIABLES HERE !!!
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

% Check input

p = inputParser;

p.addParamValue('Ewein', @(x) isscalar(x) && isstruct(x));
p.addParamValue('nekidx', [], @(x) isvector(x) && isnumeric(x));
p.addParamValue('weightfun', [], @(x) isscalar(x) && isa(x, 'function_handle'));
p.addParamValue('switching', [], @(x) isvector(x) && isnumeric(x) && all(x <= 2 & x>= 0));
p.addParamValue('kn', [], @(x) isvector(x) && isnumeric(x));
p.addParamValue('dsink', [], @(x) isvector(x) && isnumeric(x));
p.addParamValue('kremin', [], @(x) isvector(x) && isnumeric(x));
p.addParamValue('upwell', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('nutrient', [], @(x) size(x,2)==2);

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

default = p.UsingDefaults;
if ~isempty(default)
    missing = sprintf('%s, ', default{:});
    error('Missing input parameter for ewe biology: %s', missing(1:end-2));
end

In = mergestruct(p.Unmatched, p.Results);

% Calculate all non-physical parameters

Param = watercolumnecosystemsetup(In.Ewein, In.nekidx, In.weightfun, ...
        In.nutrient, In.kn, In.kremin, In.switching, Grd.z, Grd.zp);
    
% Names of functional groups

names = cell(Param.ngroup,3);
names(:,1) = cellstr(num2str((1:Param.ngroup)', 'FG%02d'));
names(:,2) = Param.names;
[names{:,3}] = deal('mmol N/m^3');

% names = cellfun(@(x) regexprep(x, '\W', ''), Param.names, 'uni', 0);
% names = reshape(names, 1, []);

% Assign parameters needed for mixed_layer

bio = Param.bio;

Biovars = rmfield(Param, {'names', 'Ewein', 'Ep', 'isnekton', 'bio', ...
                          'other'});
Biovars.loss = Param.other;

% Sinking rate (m/s, negative down) for detritus and upwelling rate for
% nutrients

Biovars.wsink = zeros(Grd.nz, Biovars.ngroup);
Biovars.wsink(:, Biovars.type == 5) = repmat(In.dsink, Grd.nz, 1);
Biovars.wsink(:,Biovars.type == 1) = In.upwell;

% Don't mix nektonic groups

ismixed = ~Param.isnekton;

% No deep water forcing

bottomval = nan(Biovars.ngroup,1);

% Diagnostics: Source and sink fluxes for each functional group

diagnostics = {'dbdt', 'ps', 'growth', 'egest', 'excrete', 'natloss', 'remin'};
if isfield(In, 'diagnostic')
    Biovars.keepdiag = ismember(diagnostics, In.diagnostic);
else
    Biovars.keepdiag = false(size(diagnostics));
    diagnames = cell(0);
end

diagtokeep = diagnostics(Biovars.keepdiag);

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
        for ig1 = 1:Biovars.ngroup
            for ig2 = 1:Biovars.ngroup
                count = count + 1;
                diagnames{count,1} = sprintf('%s%02d%02d', diagtokeep{id}, ig1, ig2);
                diagnames{count,2} = sprintf('%s, %s to %s', diagtokeep{id}, names{ig1,2}, names{ig2,2});
                diagnames{count,3} = 'mmol N m^-2 s^-1';
%                 diagnames{count} = sprintf('%s_%s_%s', diagtokeep{id}, names{ig2}, names{ig1});
            end
        end
    end
end


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Set up input parameters for ode function

WceParam = rmfield(Biovars, {'ngroup', 'wsink', 'loss', 'keepdiag'});
WceParam.z = -z;
WceParam.dz = dz;
WceParam.irr = meanqi;
WceParam.temp = temperature;
WceParam.other = Biovars.loss;
      
% Convert bio from concentration/m^3 to conc/m^2

oldbio = bsxfun(@times, oldbio, dz);
      
% Integrate using Runge-Kutta solver over the full timestep

odefun = @(t,b,P) watercolumnecosystemode(t, b, P); % specify 3 rather than 4 inputs
[tout, newbio, db, WceOut] = odewrap(@ode4, odefun, [t t+dt], oldbio, [], WceParam);           

[newbio, diag{1}] = endonly(newbio, db);
diag = [diag; struct2cell(WceOut(end))]';

% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

shouldbenan = false(size(oldbio));
shouldbenan(2:end, Biovars.type == 4) = true;

Opt = odeset('nonnegative', 1:numel(oldbio));
if any((isnan(newbio(:)) & ~shouldbenan(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout, newbio, db, WceOut] = odewrap(@ode45, odefun, [t t+dt], oldbio, Opt, WceParam);  
    [newbio, diag{1}] = endonly(newbio, db);
    diag = [diag; struct2cell(WceOut(end))]';
end
   
% Reformat diagnostics

diag = diag(Biovars.keepdiag);
for id = 1:length(diag)
    if ndims(diag{id}) == 3
        diag{id} = permute(diag{id}, [3 1 2]);
        diag{id} = reshape(diag{id}, size(diag{id},1), []);
    end
end
diag = cat(2, diag{:});

% Convert back to m^-3

newbio = bsxfun(@rdivide, newbio, dz);


%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = Biovars.wsink;