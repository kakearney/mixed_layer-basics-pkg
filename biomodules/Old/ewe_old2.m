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

%----------------------------
% Define state variables: EwE 
% functional groups plus one 
% nutrient pool
%----------------------------

names = [In.Ewein.name; {'Nutrient'}];
Biovars.ngroup = length(names);

% Strip out non-alphanumeric characters to fit fieldname requirements

names = cellfun(@(x) regexprep(x, '\W', ''), names, 'uni', 0);
names = reshape(names, 1, []);

%----------------------------
% Run Ecopath to get initial
% parameters
%----------------------------

% EwE models are parameterized on a per area basis, while the mixed layer
% model operates on a per volume basis.  All feeding calculations will
% stay in per area units, but are converted here to units more similar to
% traditional biogeochemical models. Input units are M = tons wet weight, A
% = km^2, T = year. Convert to M = mmol N, A = m^2, T = s.

Ewein = ecopathinputcheck(In.Ewein, true);
Ewein = eweunitconvert('tons wet weight/km^2/yr', 'mmol N/m^2/sec', Ewein, 'wwCfrac', .01);
% Ewein = eweunitconvert(Ewein, 'tons wet weight', 'mmol N', 'km^2', 'm^2', 'year', 'sec', 'wwCfrac', .01);

% Run Ecopath to get initial conditions

Ep = ecopathlite(Ewein);

% Tag each group as nutrient (1), producer (2), planktonic (3), nektonic
% (4), or detrital (5).  Note that the nutrient group is always the last
% group.

isnekton = false(Biovars.ngroup,1);
isnekton(In.nekidx) = true;

Biovars.type = ones(Biovars.ngroup,1);
Biovars.type(Ewein.pp == 1) = 2;
Biovars.type(Ewein.pp == 2) = 5;
Biovars.type(Ewein.pp == 0 & isnekton(1:end-1)) = 4;
Biovars.type(Ewein.pp == 0 & ~isnekton(1:end-1)) = 3;

%----------------------------
% Set initial profiles 
% (mmol N m^-2)
%----------------------------

bio = nan(Grd.nz, Biovars.ngroup);

isnek = Biovars.type == 4;
isnut = Biovars.type == 1;
isplank = ~isnek & ~isnut;

% Nektonic groups can move anywhere in the water column.  The use of row 1
% to hold their biomass value is arbitrary.  Units are mmol N m^-2.

bio(1, In.nekidx) = Ep.b(In.nekidx);

% The initial profile for the nutrient group is set by the user.  Even
% though I'll need the units to be in mmol N m^-3 in photosynthesis
% calculations, here I'll convert it to mmol N m^-2 to keep consistancy
% with the other groups.

nut = interp1(In.nutrient(:,1), In.nutrient(:,2), Grd.z);   % mmol N m^-3
dz = diff(-Grd.zp);
bio(:, isnut) = nut.*dz;                                    % mmol N m^-2

% Planktonic groups are distributed according to the weight function
% provided by the user, so that the integral over depth equals the
% ecopath-balanced value.

% TODO: right now detritus is being treated as a planktonic group. This
% should probably resemble a sediment trap profile or something like that
% though.

zmax = -(Grd.zp(end));
ztemp = -Grd.zp;

weighttot = quad(In.weightfun, 0, zmax);
weightlayer = arrayfun(@(a,b) quad(In.weightfun, a, b), ztemp(1:end-1), ztemp(2:end));
weight = weightlayer./weighttot;    % Fraction of total biomass found in each layer

biolayer = weight * Ep.b(isplank)';
bio(:,isplank) = biolayer;

biotemp = bio;
biotemp(isnan(biotemp)) = 0;
biofrac = bsxfun(@rdivide, biotemp, sum(biotemp,1));

% For the mixing routines to work properly, state variable concentrations
% must be in mmol N m^-3.  TODO: Is all this converting back and forth
% going to lead to too much numerical drift?

bio = bsxfun(@rdivide, bio, dz);

%----------------------------
% Ecosim functional response 
% parameters (these 
% parameters are not defined 
% for nutrient group)
%----------------------------

% Handling time (s)

Biovars.h = ecosimhandling(Ep.qb, Ewein.qbmaxqb0);

% Vulnerability (s^-1)

q0 = Ep.q0;
q0(:, Ewein.pp==2) = 0; % Flow to detritus is NOT consumption

Biovars.v = ecosimvulnerability(Ewein.kv, q0, Ep.b);

% Switching factor (no unit)

Biovars.sw = In.switching;

% Categorize feeding

bothnek   = Ewein.dc > 0 & bsxfun(@and, isnek(1:end-1), isnek(1:end-1)');
predisnek = Ewein.dc > 0 & bsxfun(@and, ~isnek(1:end-1), isnek(1:end-1)');
bothplank = Ewein.dc > 0 & ~bothnek & ~predisnek;

Biovars.feedtype = zeros(size(Ewein.dc));
Biovars.feedtype(bothnek)   = 1;
Biovars.feedtype(predisnek) = 2;
Biovars.feedtype(bothplank) = 3;

% Search factor (m^2 mmol^-1 s^-1) and scaling factor (no unit)

[Biovars.a, Biovars.k] = ecosimsearchlayered('type2', Ep.b, q0, Biovars.v, Biovars.h, Biovars.sw, biofrac, Biovars.feedtype);

Biovars.a(isnan(Biovars.a) | isinf(Biovars.a)) = 0;
Biovars.k(isnan(Biovars.k) | isinf(Biovars.k)) = 0;

% Function for feeding

Biovars.funcresfun = functionalresponse('type2', 'inline', 'B', 'P', 'a', 'h', 'v');

%----------------------------
% Other Ecosim-related
% parameters (also not 
% defined for nutrient)
%----------------------------

% Growth efficiency (no units)

Biovars.ge = Ep.ge;

% Fraction egested (no units)

Biovars.gs = Ewein.gs;

% Loss due to non-predatory causes (s^-1)

Biovars.loss = Ep.otherMortRate;

%----------------------------
% Producer parameters
%----------------------------

% Half-saturation constant for producers (mmol m^-3)

Biovars.kn = zeros(Biovars.ngroup,1);
Biovars.kn(Biovars.type == 2) = In.kn;

%----------------------------
% Detritus parameters
%----------------------------

% Sinking rate (m/s, negative down) for detritus and upwelling rate for
% nutrients

Biovars.wsink = zeros(Grd.nz, Biovars.ngroup);
Biovars.wsink(:, Biovars.type == 5) = repmat(In.dsink, Grd.nz, 1);
Biovars.wsink(:,Biovars.type == 1) = In.upwell;

% Remineralization rate (s^-1)

Biovars.kremin = In.kremin;

% Fraction of detritus from each group (egestion and loss) going to each
% detritus group

Biovars.df = Ewein.df;

%----------------------------
% Variables for ecosim 
% production
%----------------------------

% Biovars.rprod = Ewein.maxrelpb .* Ep.pb;
% Biovars.hprod = (Ewein.maxrelpb - 1)./Ep.b;

%----------------------------
% Mixing/forcing flags
%----------------------------

% Don't mix nektonic groups

ismixed = true(Biovars.ngroup,1);
ismixed(isnek) = false;

% No deep water forcing

bottomval = nan(Biovars.ngroup,1);

%----------------------------
% Diagnostics
%----------------------------

% Source and sink fluxes for each functional group

diagnostics = {'dbdt', 'ps', 'growth', 'excrete', 'egest', 'natloss', 'predloss', 'remin'};

count = 0;
for id = 1:length(diagnostics)
    for ig = 1:Biovars.ngroup
        count = count + 1;
        diagnames{count} = sprintf('%s_%s', diagnostics{id}, names{ig});
    end
end

% Predation flux between each predator/prey pair

count = 0;
for ipred = 1:Biovars.ngroup-1
    for iprey = 1:Biovars.ngroup-1
        count = count + 1;
        predfluxname{count} = sprintf('predflux_%s_%s', names{iprey}, names{ipred});
    end
end
diagnames = [diagnames predfluxname];

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Group input parameters and output diagnostics for code clarity

diag = cell(1,9);
params = {Biovars.type, Biovars.kn, -z, meanqi, temperature, dz, ...
          Biovars.a, Biovars.v, Biovars.h, Biovars.sw, Biovars.ge, ...
          Biovars.gs, Biovars.df, Biovars.loss, Biovars.feedtype, ...
          Biovars.k, Biovars.kremin, Biovars.funcresfun};
      
% Convert bio from concentration/m^3 to conc/m^2

oldbio = bsxfun(@times, oldbio, dz);
      
% Integrate using Runge-Kutta solver over the full timestep

[tout,newbio, diag{:}] = odewrap(@ode4, @biochange, [t t+dt], oldbio, [], params{:});           
[newbio, diag{:}] = endonly(newbio, diag{:});
   
% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

shouldbenan = false(size(oldbio));
shouldbenan(2:end, Biovars.type == 4) = true;

Opt = odeset('nonnegative', 1:numel(oldbio)); %, 'outputfcn', @odeplot);
if any((isnan(newbio(:)) & ~shouldbenan(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio, diag{:}] = odewrap(@ode45, @biochange, [t t+dt], oldbio, Opt, params{:});           
    [newbio, diag{:}] = endonly(newbio, diag{:});
end

% Concatenate diagnostics

diag = cat(2, diag{:});

% Convert back to m^-3

newbio = bsxfun(@rdivide, newbio, dz);

%----------------------------
% Main flux equation 
%----------------------------
                     
function [dbdt, ps, growth, excrete, egest, naturalLoss, ...
          predationLoss, remin, predflux] = ...
            biochange(time, bio, type, kn, z, irr, temp, dz, a, v, h, ...
                      sw, ge, gs, df, other, feedtype, k, kremin, ...
                      funcresfun)

if any(bio(:) < 0)
%     error('Populations cannot be negative');
end
     
testing = false;
if ~testing
    ps = photosyn(bio, type, kn, z, irr, temp, dz);
else
    ps = photosynsteady(bio, type, rprod, hprod);
end
                  
[predationLoss, predGain, predflux] = ingestion(bio, type, feedtype, a, k, v, h, sw, dz, funcresfun); 
                       
[growth, excrete, egest] = excreteegest(predGain, type, ge, gs, df, dz, bio);
                                    
naturalLoss = naturalmortality(bio, other, type, df);

remin = remineralization(bio, kremin, type);
                       
dbdt = ps + growth + excrete + egest - naturalLoss - predationLoss + remin;  

%----------------------------
% Photosynthesis
%----------------------------

function ps = photosyn(bio, type, kn, z, irr, temp, dz)

% Photosynthesis involves the nutrient group and all producer groups

isnut = type == 1;
isprod = type == 2;

nutrients = bio(:,isnut);
phyto = bio(:,isprod);

% Convert from mmol N m^-2 to mmol N m^-3

if isscalar(dz)
    nutrients = nutrients .* dz;
    phyto = phyto .* dz;
else
    nutrients = bsxfun(@rdivide, nutrients, dz);
    phyto = bsxfun(@rdivide, phyto, dz);
end

% Calculate photosynthesis rate (s^-1) for each producer at each depth

kn = kn(isprod)';
pstemp = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz);

% Calculate total source due to photosynthesis for all groups

ps = zeros(size(bio));
ps(:,isprod) = pstemp;  % s^-1
ps = ps .* bio;         % mmol m^-2 s^-1

% Gain in phyto is a loss in nutrient

totps = nansum(ps,2);
ps(:, type == 1) = -totps;
ps(isnan(ps)) = 0;

%----------------------------
% Photosynthesis (test)
%----------------------------

% function ps = photosynsteady(bio, type, rprod, hprod)
% 
% isprod = type == 2;
% isnut = type == 1;
% 
% bio(isnan(bio)) = 0;
% totbio = sum(bio, 1);
% 
% r = zeros(1, size(bio,2));
% h = zeros(1, size(bio,2));
% r(isprod) = rprod(isprod);
% h(isprod) = hprod(isprod);
% 
% r = repmat(r, size(bio,1), 1);
% h = repmat(h, size(bio,1), 1);
% 
% ps = r.*bio./(1 + h.*bio);
% 
% ps(:, isnut) = -sum(ps,2);

%-----------------------------
% Ingestion of prey
%-----------------------------

function [preyLoss, predGain, predflux] = ingestion(bio, type, feedtype, a, k, v, h, sw, dz, funcresfun) 

% Calculate ingestion of each prey by each predator

[ingin, ingout] = ecosimfeedlayered(funcresfun, bio(:,1:end-1), a, k, v, h, sw, feedtype);

% Separate ingestion into source for predator and sink for prey

predGain = zeros(size(bio));
preyLoss = zeros(size(bio));
predGain(:,1:end-1) = squeeze(sum(ingin, 1))';
preyLoss(:,1:end-1) = squeeze(sum(ingout,2))';

% For diagnostics, keep flux between each group (using nekton-plankton
% spread over layers)

% predflux = zeros(size(bio,2)-1, size(bio,2)-1, size(bio,1));
predflux = ingin;                                           % prey x pred x depth
predflux = permute(predflux, [3 1 2]);                      % depth x prey x pred
predflux = reshape(predflux, size(predflux,1), []);         % depth x (prey*pred)

% % Nekton-eat-nekton: straightforward use of functional response
% 
% bnn = bio(1,1:end-1)';
% 
% ingNekNek = ecosimfeed('type2', bnn, a, k, v, h, sw);
% ingNekNek(feedtype ~= 1) = 0;
% 
% % Nekton-eat-plankton: sum plankton groups over all depths first
% 
% bnp = nansum(bio(:,1:end-1), 1);                 % total over water column
% bioFrac = bsxfun(@rdivide, bio(:,1:end-1), bnp); % fraction in each layer
% 
% ingNekPlank = ecosimfeed('type2', bnp, a, k, v, h, sw);
% ingNekPlank(feedtype ~= 2) = 0;
% 
% % Plankton-eat-plankton: loop over each depth layer
% 
% for idepth = 1:size(bio,1)
%     bpp = bio(idepth,1:end-1)';
%     ingPlankPlank{idepth} = ecosimfeed('type2', bpp, a, k, v, h, sw);
%     ingPlankPlank{idepth}(feedtype ~= 3) = 0;
% end
% ingPlankPlank = cat(3, ingPlankPlank{:});
%     
% % Separate ingestion into source for predator and sink for prey
% 
% [ndepth, ng] = size(bio);
% 
% lossnn = zeros(ndepth,ng-1);        % Losses and gains both go to row 1 
% gainnn = zeros(ndepth,ng-1);
% lossnn(1,:) = sum(ingNekNek, 2)';
% gainnn(1,:) = sum(ingNekNek, 1);
% 
% losspp = sum(ingPlankPlank, 2);     % Losses and gains already per layer
% gainpp = sum(ingPlankPlank, 1);
% losspp = squeeze(losspp)';
% gainpp = squeeze(gainpp)';
% 
% gainnp = zeros(ndepth,ng-1);        % Gains stay in row 1, losses distribute over layers
% gainnp(1,:) = sum(ingNekPlank, 1);       
% lossnp = bsxfun(@times, bioFrac, sum(ingNekPlank, 2)');
% lossnp(isnan(lossnp)) = 0;
% 
% preyLoss = zeros(size(bio));
% preyLoss(:,type~=1) = lossnn + losspp + lossnp;
% predGain = zeros(size(bio));
% predGain(:,type~=1) = gainnn + gainpp + gainnp;
% 
% % For diagnostics, keep flux between each group
% 
% predflux = ingPlankPlank;
% predflux(:,:,1) = predflux(:,:,1) + ingNekNek + ingNekPlank; % prey x pred x depth
% predflux = permute(predflux, [3 1 2]);                       % depth x prey x pred
% predflux = reshape(predflux, size(predflux,1), []);          % depth x (pred*prey)

%-----------------------------
% Growth, excretion, and 
% egestion
%-----------------------------

function [growth, excrete, egest] = excreteegest(eaten, type, ge, gs, df, dz, bio);

% Remove NaN placeholders for gs (producers and detrital groups) since it
% just makes things more difficult.  

gs(isnan(gs)) = 0;

% Separate ingested food into stuff that contributes to growth, stuff
% egested (to detritus), and stuff excreted (to nutrient pool).

togrowth   = bsxfun(@times, eaten(:,type~=1), ge');
todetritus = bsxfun(@times, eaten(:,type~=1), gs');
tonutrient = bsxfun(@times, eaten(:,type~=1), (1 - ge - gs)');

% Stuff going to growth stays in the same depth and group cell

growth = zeros(size(eaten));
growth(:,type~=1) = togrowth;

% Egested stuff goes to the detritus group(s) specified by the detritus
% fate (df) parameter.  Plankton egest to their respective layers.  
% TODO: Nekton?  Should I distribute their toilet throughout the
% entire water column?  Just at the surface?  Right now I'm leaving it at
% the surface, since that's where I'm holding numbers for nekton, but it's
% a pretty arbitrary choice.   

isplank = type == 2 | type == 3;
plankbio = bio(:,isplank);
plankbio = sum(plankbio,2);
plankbiofrac = plankbio./sum(plankbio);
hasplank = plankbiofrac > .01;
deepestplank = find(hasplank, 1, 'last');
nekfrac = dz;
nekfrac(deepestplank+1:end) = 0;
nekfrac = nekfrac./(sum(nekfrac));

todetritus(:,type == 4) = nekfrac * todetritus(1, type == 4);
tonutrient(:,type == 4) = nekfrac * tonutrient(1, type == 4);

detidx = find(type == 5);

egest = zeros(size(eaten));
for igroup = 1:size(eaten,2)-1
    detfrac = zeros(1, size(eaten,2));
    detfrac(detidx) = df(igroup,:);
    egest = egest + todetritus(:,igroup) * detfrac;
end

% Excretion all goes to the nutrient pool (mmol m^-3 s^-1)

excrete = zeros(size(eaten));
excrete(:, type==1) = sum(tonutrient, 2);

%-----------------------------
% Non-predatory loss
%-----------------------------

function loss = naturalmortality(bio, lossrate, type, df)

% Calculate loss from each group

loss = zeros(size(bio));
loss(:,type~=1) = bsxfun(@times, bio(:,type~=1), lossrate');
loss(isnan(loss)) = 0;

% Loss from groups goes to detritus groups

detidx = find(type == 5);
detgain = zeros(size(bio));
for igroup = 1:size(bio,2)-1
    detfrac = zeros(1, size(bio,2));
    detfrac(detidx) = df(igroup,:);
    detgain = detgain + loss(:,igroup) * detfrac;
end

% Final loss (det groups have negative values indicating a gain)

loss = loss - detgain;

%-----------------------------
% Remineralization
%-----------------------------

function remin = remineralization(bio, kremin, type);

detremin = zeros(size(bio));
detremin = bio(:, type == 5) .* kremin;
sumremin = sum(detremin, 2);

remin(:, type == 5) = -detremin;
remin(:, type == 1) = sumremin;

%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = Biovars.wsink;