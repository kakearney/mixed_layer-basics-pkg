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
% model operates on a per volume basis.  Calculate Ecopath parameters in
% both sets of units.  Input units are M = tons wet weight, A = km^2, T =
% year. Convert to M = mmol N, A = m^2, T = s for area and M = mmol N, A =
% m^3, T = s for volume.

% First, convert the per area units to new set of per area units

fields = {...
    'ngroup'        'no unit'   
    'nlive'         'no unit'
    'ngear'         'no unit'
    'name'          'other'
    'pp'            'no unit'
    'areafrac'      'no unit'
    'bh'            'M/A'
    'pb'            '1/T'
    'qb'            '1/T'
    'ee'            'no unit'
    'ge'            'no unit'
    'gs'            'no unit'
    'dtImp'         'M/A/T'
    'b'             'M/A'
    'dc'            'no unit'
    'df'            'no unit'
    'immig'         'M/A/T'
    'emig'          'M/A/T'
    'emigRate'      '1/T'
    'ba'            'M/A/T'
    'baRate'        '1/T'
    'fleetname'     'other'
    'landing'       'M/A/T' 
    'discard'       'M/A/T'
    'discardFate'   'no unit'
    'maxrelpb'      'no unit'
    'maxrelfeed'    'no unit'
    'feedadj'       'no unit' 
    'fracsens'      'no unit'    
    'predeffect'    'no unit'
    'densecatch'    'no unit'
    'qbmaxqb0'      'no unit'
    'switchpower'   'no unit'
    'kv'            'no unit'};

fieldsByUnit = {...
    fields(strcmp(fields(:,2), 'M/A'), 1)
    fields(strcmp(fields(:,2), '1/T'), 1)
    fields(strcmp(fields(:,2), 'M/A/T'), 1)
    fields(strcmp(fields(:,2), 'no unit'), 1)
    };

ppCfrac = .05;                          % C/wet weight for phytoplankton
c2n = 17/133;                           % mol N/mol C
cmw = 12.0107;                          % molecular weight of C
ww2mm = (ppCfrac * c2n /cmw) * 1e9;     % mmol N/ton wet weight
kmsq2msq = 1e6;                         % m^2/km^2
yr2sec = 86400 * 365;                   % sec/year

factor = [ww2mm/kmsq2msq; 1/yr2sec; ww2mm/(kmsq2msq * yr2sec); 1];

Ewein(1) = In.Ewein;
for iun = 1:length(fieldsByUnit)            % Loop over unit categories
    for ifld = 1:length(fieldsByUnit{iun})  %   then each field in the cat.
        fld = fieldsByUnit{iun}{ifld};
        Ewein(1).(fld) = In.Ewein.(fld) * factor(iun);
    end
end
        
% Now calculate the per volume equivalents

area2vol = 1/In.habitatDepth;

Ewein(2) = Ewein(1);
for iun = [1 3]
    for ifld = 1:length(fieldsByUnit{iun})
        fld = fieldsByUnit{iun}{ifld};
        Ewein(2).(fld) = Ewein(1).(fld) * area2vol;
    end
end

% Run Ecopath to get initial conditions

Ep = arrayfun(@ecopathlite, Ewein);

% Tag each group as nutrient (1), producer (2), planktonic (3), nektonic
% (4), or detrital (5).  Note that the nutrient group is always the last
% group.

Biovars.type = ones(Biovars.ngroup,1);
Biovars.type(Ewein(1).pp == 1) = 2;
Biovars.type(Ewein(1).pp == 2) = 5;
Biovars.type(Ewein(1).pp == 0 & In.isnekton) = 4;
Biovars.type(Ewein(1).pp == 0 & ~In.isnekton) = 3;

%----------------------------
% Set initial profiles
%----------------------------

% Planktonic groups (including producers) start out evenly distributed in
% the habitat layers, using per volume units. Nektonic organisms remain in
% per area units.  Nutrients are set according to input profile.
% Detritus... TODO, figure out something that makes more sense, but using
% plankton method right now.
%
% Just a reminder, all columns of bio EXCEPT nektonic groups are in units
% of mmol N m^-3.  Nektonic columns are in mmol N m^-2; NaNs are just
% placeholders.

% All layers that are at least partially shallower than specificed habitat
% depth are considered habitat layers.

ishabitat = -Grd.zp(1:end-1) < In.habitatDepth; 

bio = nan(Grd.nz, Biovars.ngroup);
isnek = Biovars.type == 4;
isnut = Biovars.type == 1;
isplank = ~isnek & ~isnut;

bio(1, isnek) = Ep(1).b(isnek);
bio(:, isnut) = interp1(In.nutrient(:,1), In.nutrient(:,2), Grd.z);
nh = sum(ishabitat);
bio(ishabitat,isplank) = repmat(Ep(2).b(isplank)', nh, 1);
bio(~ishabitat,isplank) = 0;

%----------------------------
% Ecosim functional response 
% parameters
%----------------------------

% Handling time (s)

Biovars.h = 1./(Ep(1).qb .* Ewein(1).qbmaxqb0);

% Vulnerability (s^-1)

Biovars.v = bsxfun(@times, Ewein(1).kv, Ep(1).qb');

% Switching factor (no unit)

Biovars.sw = In.switching;

% Search factor (m^2 mmol^-1 s^-1 for nekton feeding, m^3 mmol^-1 s^-1 for
% all others) and scaling factors

bothnek   = Ewein(1).dc > 0 & bsxfun(@and, isnek(1:end-1), isnek(1:end-1)');
predisnek = Ewein(1).dc > 0 & bsxfun(@and, ~isnek(1:end-1), isnek(1:end-1)');
bothplank = Ewein(1).dc > 0 & ~bothnek & ~predisnek;

q0 = {Ep.q0};
q0{1}(:,Biovars.type==5) = 0; % Ecopath includes detritus flow as "consumption" of other groups
q0{2}(:,Biovars.type==5) = 0;
[aArea, kArea] = ecosimsearch('type2', Ep(1).b, q0{1}, Biovars.v, Biovars.h, Biovars.sw);
[aVol, kVol]   = ecosimsearch('type2', Ep(2).b, q0{2}, Biovars.v, Biovars.h, Biovars.sw);

Biovars.a = zeros(Biovars.ngroup-1);
Biovars.a(bothnek)   = aArea(bothnek);       % m^2/(mmol s)
Biovars.a(predisnek) = aArea(predisnek);     % m^2/(mmol s)
Biovars.a(bothplank) = aVol(bothplank);      % m^3/(mmol s)

Biovars.k = zeros(Biovars.ngroup-1);
Biovars.k(bothnek)   = kArea(bothnek);       % m^2/(mmol s)
Biovars.k(predisnek) = kArea(predisnek);     % m^2/(mmol s)
Biovars.k(bothplank) = kVol(bothplank);      % m^3/(mmol s)

Biovars.feedtype = zeros(size(Ewein(1).dc));
Biovars.feedtype(bothnek)   = 1;
Biovars.feedtype(predisnek) = 2;
Biovars.feedtype(bothplank) = 3;

%----------------------------
% Other Ecosim-related
% parameters
%----------------------------

% Growth efficiency (no units)

Biovars.ge = Ep(1).ge;

% Fraction egested (no units)

Biovars.gs = Ewein(1).gs;

% Loss due to non-predatory causes (s^-1)

Biovars.loss = Ep(1).otherMortRate;

%----------------------------
% Producer parameters
%----------------------------

% Half-saturation constant for producers (mmol m^-3)

Biovars.kn = zeros(Biovars.ngroup,1);
Biovars.kn(Biovars.type == 2) = In.kn;

%----------------------------
% Detritus parameters
%----------------------------

% Sinking rate (m/s, negative down) for detritus

Biovars.wsink = zeros(Grd.nz, Biovars.ngroup);
Biovars.wsink(:, Biovars.type == 5) = repmat(In.dsink, Grd.nz, 1);

% Remineralization rate (s^-1)

Biovars.kremin = In.kremin;

% Fraction of detritus from each group (egestion and loss) going to each
% detritus group

Biovars.df = Ewein(1).df;

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

diagnostics = {'dbdt', 'ps', 'growth', 'excrete', 'egest', 'natloss', 'predloss'};

% All diagnostics are returned for each group in the model

count = 0;
for id = 1:length(diagnostics)
    for ig = 1:Biovars.ngroup
        count = count + 1;
        diagnames{count} = sprintf('%s_%s', diagnostics{id}, names{ig});
    end
end

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Group input parameters and output diagnostics for code clarity

diag = cell(1,7);
params = {Biovars.type, Biovars.kn, -z, meanqi, temperature, dz, ...
          Biovars.a, Biovars.v, Biovars.h, Biovars.sw, Biovars.ge, ...
          Biovars.gs, Biovars.df, Biovars.loss, Biovars.feedtype, Biovars.k};
      
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


% [tout,newbio,diag{:}] = odewrap(@ode45, @biochange, [t t+dt], oldbio, Opt, ...
%                                  params{:});
% [newbio, diag{:}] = endonly(newbio, diag{:});

% Concatenate diagnostics

diag = cat(2, diag{:});

%----------------------------
% Main flux equation 
%----------------------------
                     
function [dbdt, ps, growth, excrete, egest, naturalLoss, ...
          predationLoss] = biochange(time, bio, type, kn, z, irr, temp, ...
                                     dz, a, v, h, sw, ge, gs, df, other, ...
                                     feedtype, k)

if any(bio(:) < 0)
%     error('Populations cannot be negative');
end
                                 
ps = photosyn(bio, type, kn, z, irr, temp, dz);
                  
[predationLoss, predGain] = ingestion(bio, type, feedtype, a, k, v, h, sw, dz); 
                       
[growth, excrete, egest] = excreteegest(predGain, type, ge, gs, df, dz);
                                    
naturalLoss = naturalmortality(bio, other, type, dz);
                       
dbdt = ps + growth + excrete + egest - naturalLoss - predationLoss;  

%----------------------------
% Photosynthesis
%----------------------------

function ps = photosyn(bio, type, kn, z, irr, temp, dz)

% Calculate photosynthesis rate (s^-1) for each producer at each depth

isnut = type == 1;
isprod = type == 2;

nutrients = bio(:,isnut);
phyto = bio(:,isprod);
kn = kn(isprod)';
pstemp = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz);

% Calculate total source due to photosynthesis for all groups

ps = zeros(size(bio));
ps(:,isprod) = pstemp;  % s^-1
ps = ps .* bio;         % mmol m^-3 s^-1

%-----------------------------
% Ingestion of prey
%-----------------------------

function [preyLoss, predGain] = ingestion(bio, type, feedtype, a, k, v, h, sw, dz) 

% Nekton-eat-nekton: biomass is in  mmol m^-2 and a is in m^2/(mmol s),
% ingestion is in mmol m^-2 s^-1

bnn = bio(1,1:end-1)';
hnn = h;
knn = k;
ann = a;
snn = sw;
vnn = v;

ingNekNek = ecosimfeed('type2', bnn, ann, knn, vnn, hnn, snn);
% ingNekNek = estype2(ann, bnn, vnn, hnn, knn, pnn);
ingNekNek(feedtype ~= 1) = 0;

% Nekton-eat-plankton: biomass for nekton is in mmol m^-2 but plankton is
% in mmol m^-3, a is in m^2/(mmol s).  Integrate plankoton over water
% column before calculation, so ingestion is in mmol m^-2 s^-1

bioIntegrated = bsxfun(@times, bio, dz);
planktonIntegrated = bioIntegrated(:,type~=4 & type~=1);
bioIntegrated = bioIntegrated(:,type~=1);
bioFrac = bsxfun(@rdivide, bioIntegrated, nansum(bioIntegrated));

bnp = zeros(size(a,1),1);
bnp(type == 4) = bio(1,type == 4);
bnp(type ~=4 & type ~=1) = sum(planktonIntegrated, 1);

hnp = h;
knp = k;
anp = a;
vnp = v;
snp = sw;

ingNekPlank = ecosimfeed('type2', bnp, anp, knp, vnp, hnp, snp);
% ingNekPlank = estype2(anp, bnp, vnp, hnp, knp, pnp);
ingNekPlank(feedtype ~= 2) = 0;

% Plankton-eat-plankton: biomass in mmol m^-3, a is in m^3/(mmol s),
% ingestion is in mmol m^-3 s^-1.  Do calculation for each layer.

hpp = h;
kpp = k;
app = a;
vpp = v;
spp  = sw;

for idepth = 1:size(bio,1)
    bpp = bio(idepth,1:end-1)';
    ingPlankPlank{idepth} = ecosimfeed('type2', bpp, app, kpp, vpp, hpp, spp);
%     ingPlankPlank{idepth} = estype2(app, bpp, vpp, hpp, kpp, ppp);
    ingPlankPlank{idepth}(feedtype ~= 3) = 0;
end
ingPlankPlank = cat(3, ingPlankPlank{:});
    
% Separate ingestion into source for predator and sink for prey

lossnn = sum(ingNekNek, 2)';    % 1 per group, mmol m^-2 s^-1
gainnn = sum(ingNekNek, 1);

losspp = sum(ingPlankPlank, 2); % 1 per group and layer, mmol m^-3 s^-1
gainpp = sum(ingPlankPlank, 1);
losspp = squeeze(losspp)';
gainpp = squeeze(gainpp)';

gainnp = sum(ingNekPlank, 1);               % 1 per group, mmol m^-2 s^-1
lossnp = sum(ingNekPlank, 2);               % 1 per group, mmol m^-2 s^-1
lossnp = bsxfun(@times, lossnp', bioFrac);  % 1 per group and layer, mmol m^-2 s^-1
lossnp = bsxfun(@rdivide, lossnp, dz);      % 1 per group and layer, mmol m^-3 s^-1

preyLoss = zeros(size(bio));
preyLoss(:,type~=1) = losspp + lossnp;
preyLoss(1,type~=1) = preyLoss(1,type~=1) + lossnn;
predGain = zeros(size(bio));
predGain(:, type~=1) = gainpp;
predGain(1,type~=1) = predGain(1,type~=1) + gainnn + gainnp;

%-----------------------------
% Growth, excretion, and 
% egestion
%-----------------------------

function [growth, excrete, egest] = excreteegest(eaten, type, ge, gs, df, dz);

% Remove NaN placeholders for gs (producers and detrital groups) since it
% just makes things more difficult.  

gs(isnan(gs)) = 0;

% Separate ingested food into stuff that contributes to growth, stuff
% egested (to detritus), and stuff excreted (to nutrient pool).

togrowth   = bsxfun(@times, eaten(:,type~=1), ge');
todetritus = bsxfun(@times, eaten(:,type~=1), gs');
tonutrient = bsxfun(@times, eaten(:,type~=1), (1 - ge - gs)');

% Stuff going to growth stays in the same units as feeding (a mix of mmol
% m^-3 s^-1 and mmol m^-2 s^-1, depending on the type of feeder).

growth = zeros(size(eaten));
growth(:,type~=1) = togrowth;

% Stuff going to nutrient and detritus pools need to be in mmol m^-3 s^-1,
% so need to convert for stuff coming from a nektonic group.
% TODO: But where?  Should I distribute their toilet throughout the entire
% habitat area?  Just at the surface?  Right now I'm leaving it at the
% surface, since that's where I'm holding numbers for nekton, but it's a
% pretty arbitrary choice.

todetritus(1,type==4) = todetritus(1,type==4)./dz(1);
tonutrient(1,type==4) = tonutrient(1,type==4)./dz(1);

% Excretion all goes to the nutrient pool (mmol m^-3 s^-1)

excrete = zeros(size(eaten));
excrete(:, type==1) = sum(tonutrient, 2);

% Egestion is divided between detrital groups according to the detritus
% fate parameter (mmol m^-3 s^-1)

df = permute(df, [3 1 2]);
egestfate = bsxfun(@times, todetritus, df); % prey x pred x det

egest = zeros(size(eaten));
egest(:,type == 5) = squeeze(sum(egestfate,2));

%-----------------------------
% Non-predatory loss
%-----------------------------

function loss = naturalmortality(bio, lossrate, type, dz)

% TODO This loss should probably go to a detritus group, but for now
% non-predatory loss is being removed from the system.

loss = zeros(size(bio));
loss(:,type~=1) = bsxfun(@times, bio(:,type~=1), lossrate');


%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% MODIFY CODE HERE

wsink = Biovars.wsink;