function varargout = np(action, varargin)
%NP Nutrient-phytoplankton biological module
%
% See biomodule.m for full syntax details.
%
% This module adds a simple nutrient-phytoplankton model to the
% mixed_layer model.  The model is based on the NP model described in
% Sarmiento and Gruber (2008) Chapter 4, where
%
%   dP/dt = P * (Vmax * N/(Kn + N) - lambdap)
%   dN/dt = P * (-Vmax * N/(Kn + N) + mup * lambdap)
%
% The maximum growth rate Vmax is determined by temperature and light
% limitation, following Eppley (1972) and Platt and Jasby (1976),
% respectively.  The half-saturation constant (Kn), loss rate (lambdap),
% and fraction of remineralization (mup) are provided by the user.
%
% Sarmiento J, Gruber N (2006) Ocean Biogeochemical Dynamics. Princeton
% University Press 
%
% User-specified input variables (passed to mixed_layer as parameter/value
% pairs)
%
%   n:      n x 2 depth profile of initial nutrients, where column 1 gives
%           the depth values (negative down) and column 2 holds the
%           concentrations of nutrients (mmol N m^-3)
%
%   p:      n x 2 depth profile of phytoplankton, where column 1 gives the
%           depth values (negative down) and column 2 holds the
%           concentrations of phytoplankton (mmol N m^-3)
%
%   kn:     half-saturation for nutrient uptake by phytoplankton (mmol N
%           m^-3) 
%
%   loss:   loss rate for phytoplankton (s^-1)
%
%   remin:  fraction of phytoplankton loss that is remineralized

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

p.addParamValue('n', [], @(x) size(x,2)==2);
p.addParamValue('p', [], @(x) size(x,2)==2);
p.addParamValue('kn', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('loss', [], @(x) isscalar(x) && isnumeric(x));
p.addParamValue('remin', [], @(x) isscalar(x) && isnumeric(x) && x <= 1 && x>= 0);

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

default = p.UsingDefaults;
if ~isempty(default)
    missing = sprintf('%s, ', default{:});
    error('Missing input parameter for np biology: %s', missing(1:end-2));
end

In = mergestruct(p.Unmatched, p.Results);

% Two state variables, N and P, both mixed

bio(:,1) = interp1(In.n(:,1), In.n(:,2), Grd.z);
bio(:,2) = interp1(In.p(:,1), In.p(:,2), Grd.z);

ismixed = true(1,2);

% Force deep-water nutrients, not phytoplankton

bottomval = [In.n(end,2) NaN];  

% Extra variables

Biovars.kn = In.kn;
Biovars.loss = In.loss;
Biovars.remin = In.remin;

% Names

names = {...
   'N', 'nutrient',         'mmol N m^-3'
   'P', 'phytoplankton',    'mmol N m^-3'};

% Diagnostics

diagnames = {...
    'ps',       'photosynthesis',   'mmol N m^-3 s^-1'
    'ploss'     'phyto mortality',  'mmol N m^-3 s^-1'
    'premin',   'remineralization', 'mmol N m^-3 s^-1'};

% diagnames = {'psPerBio', 'rate of photosynthesis per unit P', 's^-1'};

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);

                       
% Integrate using Runge-Kutta solver over the full timestep

[tout,newbio,dbdt,ps,ploss,premin] = ...
                odewrap(@ode4, @biochange, [t t+dt], oldbio, [], ...
                        temperature, meanqi, -z, dz, Biovars.remin, ...
                        Biovars.loss, Biovars.kn);

[newbio,dbdt,ps,ploss,premin] = endonly(newbio,dbdt,ps,ploss,premin);

% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio,dbdt,ps,ploss,premin] = odewrap(@ode45, @biochange, [t t+dt], oldbio, [], ...
                        temperature, meanqi, -z, dz, Biovars.remin, ...
                        Biovars.loss, Biovars.kn);
    [newbio,dbdt,ps,ploss,premin] = endonly(newbio,dbdt,ps,ploss,premin);
end

% If still no good, error

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end

diag = [ps ploss premin];
               
%-----------------------------
% Rate of change function
%-----------------------------

function [dbdt, ps, ploss, premin] = biochange(t, bio, temp, irr, z, dz, remin, loss, kn)

nutrients = bio(:,1);
phyto     = bio(:,2);

psRate = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz);
ps = psRate .* phyto;
ploss = loss .* phyto;
premin = remin .* loss .* phyto;

% Uptake of nutrient due to photosynthesis

dbdt = zeros(size(bio));

dbdt(:,1) = premin - ps; % phyto .* (remin .* loss - psRate);
dbdt(:,2) = ps - ploss; % phyto .* (psRate - loss);


%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = zeros(size(oldbio));



