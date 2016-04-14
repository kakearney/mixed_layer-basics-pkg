function varargout = wce(action, varargin)
%WCE Water column ecosystem biological module
%
% See biomodule.m for full syntax details.
%
% This module simulates a mixed planktonic-nektonic ecosystem, based on
% a combination of models derived mainly from NEMURO and Kerim Aydin's
% version of Ecosim, with a little bit of COBALT thrown in for flavor.
%
% A note on units: The biomass of all critters is saved to file in mol
% N[Si][Fe]/m^3.  For nektonic critters, all biomass is placed in the
% surface cell, and actually represents the total over the entire water
% column; multiply by the thickness of the surface layer to get the true
% biomass, in mol N/m^2.
%
% User-specified input variables (passed to mixed_layer as parameter/value
% pairs):
%
%   Ewein:      1 x 1 Ewe input structure.  Must include all ecopath fields
%               (see ecopathlite) in units of mol N/m^2/s, as well as the
%               following fields:
%                       
%               x:      Vulnerability parameter.  Can be either a ngroup x
%                       1 vector of log-transformed anomaly-from-base
%                       values (same as P4/P5 inputs in aydin-ecosim, range
%                       from -Inf to Inf w/ default of 0, xi in equation
%                       below),  or an ngroup x ngroup array of
%                       non-tranformed values for each predator-prey pair
%                       (range from 1 to Inf w/ default of 2, Xij in
%                       equation below). 
%
%                       Xij = exp(xi + xj) + 1                               
%
%               d:      Handling time parameter, same format as x. 
%
%               theta:  Switching parameter.  Same format as x but with
%                       slightly different transform  
%
%                       THij = exp(0.5 * (thi + thj))
%
%                       TH = 1 yields a type 2 functional response
%                       TH = 2 yields a type 3 functional response 
%
%   NemParam:   1 x 1 structure of NEMURO parameters (see
%               nemuroinputparser) 
%
%   types:      ngroup x 1 cell array of strings, indicating what type of
%               critter each functional group is. 
%               'n':    nekton, not affected by physical mixing, feed over
%                       entire water column 
%               'z':    zooplankton, mixed, feed only in the layer where
%                       they are located 
%               Can also be any of the 11 nemuro state variables ('ps',
%               'pl', 'zs', 'zl', 'zp', 'no3', 'nh4', 'pon', 'don',
%               'sioh4', 'opal'), all of which are planktonic.
%
%   bnem0:      n x 12 array.  Column 1 holds depth values (m, negative
%               down) and the remaining columns hold the initial biomass
%               values for nemuro-derived variables (mol/l).  Currently
%               only the non-living profiles are used, but the entire
%               matrix is used as input for consistency with the NEMURO
%               module.
%
%   nomix:      logical scalar, if true bioligcal variables are not subject
%               to mixing (for debugging purposes) [false] 
%
%   ecosimpp:   logical scalar, if true ecosim primary production function
%               is used, decoupling biology from nutrient constraints (for
%               debuging purposes) [false]
%
%   mld:        mixed layer depth, used for initial plankton distributions
%               and per-area-to-per-volume rate conversions (m, negative
%               down) [-50]
%
%   temp:       temperature associated with initial mass-balanced values,
%               used for grazing functional response (deg C) [0] 
%
%   kgra:       nexz x 1 array, temperature coefficient associated with
%               non-nemuro zooplankton groups (i.e. those with type 'z').
%               Default is 0.0693 deg C^-1 for all, i.e. a Q10 of 2.
%
%   Kfe:        1 x 2 array, half-saturation constants for iron uptake for
%               small and large phytoplankton, respectively (molFe/m^3)
%               [6e-7 3e-6]
%
%   kfe2n:      1 x 2 array, half-saturation constants for internal Fe:N
%               ratio for small and large phytoplankton, respectively
%               (molFe/molN) [[6.625e-05 0.0001325], i.e. [10 20]
%               umolFe/molC assuming Redfield]
%
%   fe2nmax:    1 x 2 array, maximum internal Fe:N ratio for small and
%               large phytoplankton, respectively (molFe/molN) 
%               [[0.00033125 0.0033125], i.e. [50 500] umolFe/molC assuming
%               Redfield]
%
%   fe2nupfac:  scalar, Fe:N uptake ratio (molFe/molN) [100e-6]
%
%   ligbkg:     Ligand background concentration (mol/m^3) [1.0e-6]
%
%   alphascav:  Iron scavenging coefficient (s^-1) [1.5855e-06, i.e. 50
%               yr^-1]
%
%   remineff:   Fraction of particulate iron remineralization relative to
%               organic (N) material (no unit) [0.25]
%
%   kiscav:     Half-saturation constant for light's influence on ligand
%               binding (W/m^2) [1.0]
%
%   kliglo:     Lower limit of ligand binding under low-light conditions
%               (m^3/mol) [3.0e8]
%
%   klighi:     Upper limit of ligand binding under high-light conditions
%               (m^3/mol) [1.0e5]  
%
%   m0exp:      Exponent for mortality function, of form M0 = aB^(m0exp). A
%               value of 1 leads to linear mortality and 2 to quadratic
%               mortality. Can also be a nlive x 1 array of values to allow
%               different functions for each critter.[2] 
%
%   reroute:    n x 5 cell array.  This allows you to reroute fluxes from
%               the original path defined in the nemuro model.  Each row
%               desribes as change, with column as follows:
%               col 1:  name of flux (gpp, gra, pre, res, exx, exc, ege, 
%                       mor, dec)
%               col 2:  name of original source group
%               col 3:  name of original sink group
%               col 4:  name of new sink group
%               col 5:  fraction of flux to reroute  
%
%   odesolver:  cell array of strings, indicating which solvers to use.  If
%               the first one fails to integrate a timestep (i.e. causes
%               something to become negative, NaN, or Inf), the next one is
%               tried.  See integratebio.m for choices. [{'euler'}]
%
%   diapause:   logical scalar, if true, a portion of the ZL group will
%               migrate to deep water seasonally, according to the
%               parameters below
%
%   dday:       1 x 3 array indicating day of year to 1) begin resurfacing,
%               2) stop resurfacing and recombine diapause/no-diapause
%               groups, and 3) migrate to deep water
%
%   dfrac:      fraction of ZL population that seasonally migrates
%
%   eatdcop:    logical scalar.  If true, predators continue to feed on
%               migrates copepods.  If false, the diapausing fraction of
%               the ZL group is ignored by flux calculations.


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

%------------------------------
% Parse and check input
%------------------------------

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParamValue('nomix', false, @(x) islogical(x) && isscalar(x));
p.addParamValue('ecosimpp', false, @(x) islogical(x) && isscalar(x));
p.addParamValue('Ewein', [], @isstruct);
p.addParamValue('NemParam', [], @isstruct);
p.addParamValue('types', []);
p.addParamValue('bnem0', []);
p.addParamValue('odesolver', {'euler'});
p.addParamValue('mld', -50, @(x) isscalar(x) && x<=0);
p.addParamValue('temp', 0, @(x) isscalar(x));
p.addParamValue('kgra', 0.0693, @(x) isvector(x));
p.addParamValue('reroute', cell(0,5), @(x) isempty(x) || (iscell(x) && size(x,2)==5));
% p.addParamValue('Fe', [], @isstruct);
% p.addParamValue('RCN',       6.625,               @(x) isnumeric(x) && iscalar(x));
% p.addParamValue('RCFe',      33e3,                @(x) isnumeric(x) && iscalar(x));
% p.addParamValue('Kfe',       [0.033e-6 0.188e-6], @(x) isnumeric(x) && isvector(x) && length(x)==2);
% p.addParamValue('Kfe',       [6e-7 3e-6], @(x) isnumeric(x) && isvector(x) && length(x)==2);
% p.addParamValue('kfe2n',     [10.0 20.0]*1e-6*106/16,   @(x) isnumeric(x) && isvector(x) && length(x)==2);
% p.addParamValue('fe2nmax',   [50.0 500.0]*1e-6*106/16,  @(x) isnumeric(x) && isvector(x) && length(x)==2);
% p.addParamValue('fe2nupfac', 100e-6,                    @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('ligbkg',    1.0e-6,                    @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('alphascav', 50.0/(86400*365),          @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('remineff',  0.25,                      @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('kiscav',    1.0,                       @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('kliglo',    3e8,                       @(x) isnumeric(x) && isscalar(x));
% p.addParamValue('klighi',    1e5,                       @(x) isnumeric(x) && isscalar(x));
p.addParamValue('remineff',  0.5,                       @(x) isnumeric(x) && isscalar(x));
p.addParamValue('kliglo',    1.0e12,                    @(x) isnumeric(x) && isscalar(x));
p.addParamValue('klighi',    1.0e8,                     @(x) isnumeric(x) && isscalar(x));
p.addParamValue('iofescav',  10,                        @(x) isnumeric(x) && isscalar(x));
p.addParamValue('alphascav', 15.0,                      @(x) isnumeric(x) && isscalar(x));
p.addParamValue('felig2don', 0,                         @(x) isnumeric(x) && isscalar(x));
p.addParamValue('ligbkg',    1e-9,                      @(x) isnumeric(x) && isscalar(x));
p.addParamValue('kfe2n',     [3.0 6.0]*1e-6*106/16,     @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('Kfe',       [1e-10 5e-10],             @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('fe2nmax',   [50.0 500.0]*1e-6*106/16,  @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('fe2nupfac', 15e-6,                     @(x) isnumeric(x) && isscalar(x));
p.addParamValue('m0exp',     2,                         @(x) isnumeric(x) && (isscalar(x) || isvector(x)));
p.addParamValue('diapause',  false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('dday',      [120 134 243],             @(x) isnumeric(x) && isvector(x) && length(x)==3);
p.addParamValue('dfrac',     0.9,                       @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParamValue('eatdcop',   true,                      @(x) islogical(x) && isscalar(x));

p.parse(In);
BioIn = p.Results;

nodefault = {'Ewein', 'NemParam', 'types', 'bnem0'};
tf = ismember(p.UsingDefaults, nodefault);
missing = sprintf('%s, ', p.UsingDefaults{tf});
if any(tf)
    error('Missing input for wce module: %s', missing(1:end-2));
end


%------------------------------
% Variable names
%------------------------------

  % State variables
  % ---------------

% Variables with 1:1 correspondance to nemurokak.  All of these will appear
% in the final model.

nemnames = {
    'PS',       'Small Phytoplankton',          'mol N m^-3' 
    'PL',       'Large Phytoplankton',          'mol N m^-3'
    'ZS',       'Small Zooplankton',            'mol N m^-3'
    'ZL',       'Large Zooplankton',            'mol N m^-3'
    'ZP',       'Predatory Zooplankton',        'mol N m^-3'
    'NO3',      'Nitrate',                      'mol N m^-3'
    'NH4',      'Ammmonium',                    'mol N m^-3'
    'PON',      'Particulate Organic Nitrogen', 'mol N m^-3'
    'DON',      'Dissolved Organic Nitrogen',   'mol N m^-3'
    'SiOH4',    'Silicate',                     'mol Si m^-3'
    'Opal',     'Particulate Opal',             'mol Si m^-3'
    'Fe'        'Dissolved Iron'                'mol Fe m^3'
    'PLsi'      'Large phytoplantkon Silica',   'mol Si m^-3'
    'PSfe'      'Small Phytoplankton Iron'      'mol Fe m^-3'
    'PLfe'      'Large Phytoplankton Iron'      'mol Fe m^-3'
    'POFe'      'Particulate iron'              'mol Fe m^-3'};

% Figure out which of the Ecopath model groups correspond to these
% nemuro-derived variables.

[tf, loc] = ismember(lower(nemnames(:,1)), BioIn.types); 

nbsv = BioIn.Ewein.ngroup + sum(~tf);

names = cell(nbsv,3);
names(loc(tf),:) = nemnames(tf,:);
names(loc(tf),2) = BioIn.Ewein.name(loc(tf));
names(BioIn.Ewein.ngroup+1:end,:) = nemnames(~tf,:);

% Name the remaining nekton and zooplankton groups, using N# and Z# for
% short names and Ecopath names for long names

isemp = cellfun('isempty', names(:,2));
names(isemp,2) = BioIn.Ewein.name(isemp);

[names{isemp,3}] = deal('mol N m^-3');
isn = strcmp(BioIn.types, 'n');
isz = strcmp(BioIn.types, 'z');

names(isn,1) = cellstr(num2str((1:sum(isn))', 'N%02d'));
names(isz,1) = cellstr(num2str((1:sum(isz))', 'Z%02d'));

  % Identifiers
  % ---------------

% Marker for nemuro-derived variables, zoo, and nekton

isnem = ismember(names(:,1), nemnames(:,1));
[blah, nemidx] = ismember(nemnames(:,1), names(:,1));

isextrazoo = regexpfound(names(:,1), 'Z[0-9][0-9]');
isnek = regexpfound(names(:,1), 'N[0-9][0-9]');
iszoo = ismember(names(:,1), {'ZS','ZL','ZP'}) | isextrazoo;
isphy = ismember(names(:,1), {'PS','PL'});

% Links and identifiers

Biovars.links = zeros(nbsv);
Biovars.links(bsxfun(@and, iszoo, iszoo')) = 1; % Z-eat-Z
Biovars.links(bsxfun(@and, iszoo, isnek')) = 2; % N-eat-Z
Biovars.links(bsxfun(@and, isnek, isnek')) = 3; % N-eat-N
Biovars.links(bsxfun(@and, isphy, iszoo')) = 4; % Z-eat-P
[r,c] = find(BioIn.Ewein.dc == 0);
idx = sub2ind(size(Biovars.links), r, c);
Biovars.links(idx) = 0;

Biovars.idx.ps    = find(strcmp(names(:,1), 'PS'));
Biovars.idx.pl    = find(strcmp(names(:,1), 'PL'));
Biovars.idx.no3   = find(strcmp(names(:,1), 'NO3'));
Biovars.idx.nh4   = find(strcmp(names(:,1), 'NH4'));
Biovars.idx.sioh4 = find(strcmp(names(:,1), 'SiOH4'));
Biovars.idx.don   = find(strcmp(names(:,1), 'DON'));
Biovars.idx.pon   = find(strcmp(names(:,1), 'PON'));
Biovars.idx.opal  = find(strcmp(names(:,1), 'Opal'));
Biovars.idx.fe    = find(strcmp(names(:,1), 'Fe'));
Biovars.idx.plsi  = find(strcmp(names(:,1), 'PLsi'));
Biovars.idx.mys   = nbsv + 1; % "Mystery" box, for source/sink terms coming from or going to nowhere
Biovars.idx.psfe  = find(strcmp(names(:,1), 'PSfe'));
Biovars.idx.plfe  = find(strcmp(names(:,1), 'PLfe'));
Biovars.idx.pofe  = find(strcmp(names(:,1), 'POFe'));
Biovars.idx.zl    = find(strcmp(names(:,1), 'ZL'));

Biovars.isnek = isnek;

Biovars.diapause = BioIn.diapause;
if Biovars.diapause
    names = [
        names
        {...
        'ZL1',       'Large Zooplankton (no diapause)', 'mol N m^-3'
        'ZL2',       'Large Zooplankton (diapause)',    'mol N m^-3'}];
    Biovars.idx.zl1 = nbsv+1;
    Biovars.idx.zl2 = nbsv+2;
    
    nbsv = nbsv + 2;
    
    linktmp = zeros(nbsv);
    linktmp(1:nbsv-2,1:nbsv-2) = Biovars.links;
    Biovars.links = linktmp;
    
    Biovars.isnek = [Biovars.isnek; false; false];
end

  % Diagnostics
  % ---------------
  
diagnames = {
    'PSlightlim',   'Light limitation (small)',         'no units'
    'PLlightlim',   'Light limitation (large)',         'no units'
    'PSno3lim',     'Nitrate limitation (small)',       'no units'
    'PLno3lim',     'Nitrate limitation (large)',       'no units'
    'PSnh4lim',     'Ammonium limitation (small)',      'no units'
    'PLnh4lim',     'Ammonium limitation (large)',      'no units'
    'PSpsmax',      'Max photosynthesis @temp (small)'  'no units'
    'PLpsmax',      'Max photosynthesis @temp (large)'  'no units'
    'PLsilim',      'Silica limitation (large)'         'no units'
    'I',            'Irradiance'                        'W m^-2'
    'kappa',        'Attenuation coefficient'           'm^-1'
    'kp'            'Attenuation self-shading only',    'm^-1'
    'PSfelim'       'Iron limitation (small)'           'no units'
    'PLfelim'       'Iron limitation (large)'           'no units'
    'PSfe2n'        'Fe:N ratio (small)'                'molFe/molN'
    'PLfe2n'        'Fe:N ratio (large)'                'molFe/molN'
    'PSfedef'       'Fe deficiency (small)'             'molFe/molN'
    'PLfedef'       'Fe deficiency (large)'             'molFe/molN'
    'extrasi'       'Extra silica due to negatives'     'mol Si m^-3'};
  
  % Intermediate fluxes
  % --------------------
  
% db/dt for each group

dbnames = names;
for idb = 1:size(dbnames,1)
    dbnames{idb,1} = ['d' dbnames{idb,1}];
    dbnames{idb,2} = ['dB/dt: ' dbnames{idb,2}];
    dbnames{idb,3} = 'molN m^-3 s^-1';
end

% Fluxes between groups, by process type

Biovars.ecosimppflag = BioIn.ecosimpp;

fluxlist = listfluxes('wce', Biovars.idx, Biovars.links);
if Biovars.ecosimppflag
    ismissing = strcmp(fluxlist(:,1), 'gpp') | ...
                strcmp(fluxlist(:,1), 'exx') | ...
                strcmp(fluxlist(:,1), 'res');
else
    ismissing = strcmp(fluxlist(:,1), 'npp');
end
fluxlist = fluxlist(~ismissing,:);

% Add user-rerouted fluxes

if isempty(BioIn.reroute)
    Biovars.reroute = cell(0,5);
else
    Biovars.reroute = BioIn.reroute;
    [tf, loc] = ismember(Biovars.reroute(:,2:4), names(:,1));
    if ~all(tf)
        error('Unrecognized critter in reroute table');
    end
    Biovars.reroute(:,2:4) = num2cell(loc);
end

fluxlist = [fluxlist; Biovars.reroute(:,[1 2 4])];
nd = size(fluxlist,1);

fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', fluxlist{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', fluxlist{id,:});
    fluxnames{id,3} = 'mol m^-3 s^-1';
end

Biovars.wceflux = fluxlist;

% Combine intermediate fluxes and extra diagnostics

diagnames = [dbnames; fluxnames; diagnames];

%----------------------------
% Initial biomass for
% state variables
%----------------------------

nz = length(Grd.z);

bio = zeros(nz, nbsv);

% Nemuro-derived non-live variables start with user-defined values in
% mol/m^3 (includes NO3 through Fe)

if all(BioIn.bnem0 >= 0)
    BioIn.bnem0(:,1) = -BioIn.bnem0(:,1); % Tired of trying to remember +up or down
end

bio(:,nemidx(6:12)) = interp1(BioIn.bnem0(:,1), BioIn.bnem0(:,7:13), Grd.z); % mol/m^3

% Plankton groups distribute their biomass evenly throughout the mixed
% layer (mld as defined by the user, not T and S)

isabovemld = Grd.z >= BioIn.mld;
nlayer = sum(isabovemld);
mld = -Grd.zp(nlayer+1);

Ewein = ecopathinputcheck(BioIn.Ewein, true);
Ep = ecopathlite(Ewein);

isplank = iszoo | isphy;

bplank = Ep.b(isplank)'./mld;   % mol N/m^3
bio(isabovemld,isplank) = repmat(bplank, nlayer, 1); 

% Nekton biomass is stored in top layer for convenience.  Although actually
% per area, here stored as per volume for consistency

bnek = Ep.b(isnek)'./-Grd.zp(2); % mol N/m^3
bio(1,isnek) = bnek;

% PL silica is proportional to PL N

pln = bio(:, Biovars.idx.pl);
plsi = pln .* BioIn.NemParam.RSiN;
bio(:, Biovars.idx.plsi) = plsi;

% PS and PL iron starts at 0

bio(:,[Biovars.idx.psfe Biovars.idx.plfe]) = 0;

% So does POFe

bio(:, Biovars.idx.pofe) = 0;

% If diapause, split ZL biomass

if Biovars.diapause
    bio(:, Biovars.idx.zl1) = bio(:, Biovars.idx.zl);
    bio(:, Biovars.idx.zl) = 0;
end
    
%----------------------------
% Indicators for mixing and 
% bottom-forcing
%----------------------------

% Plankton mixed, nekton not
% If no-mix (debugging), nothing is mixed.

if BioIn.nomix
    ismixed = false(nbsv,1);
else
    ismixed = ~isnek;
end

if Biovars.diapause
    ismixed([Biovars.idx.zl1 Biovars.idx.zl2]) = true;
end

% No bottom forcing

bottomval = nan(nbsv,1);

%----------------------------
% Variables needed for ODE
%----------------------------

% ODE solver

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

% Rerouting of fluxes (set above in flux)

% Biovars.reroute = BioIn.reroute;

% Extended NEMURO parameters with 1:1 correspondence to nemurokak
% (only defined for original-11 NEMURO-derived variables, so reassign to
% proper indices) 

Np = nemuroflexinput(BioIn.NemParam);

Biovars.alpha1    = Np.alpha1;              % m^-1  
Biovars.alpha2    = Np.alpha2/1000;         % (m^3/molN)/m
Biovars.usesteele = Np.usesteele;           % logical
Biovars.RSiN      = Np.RSiN;                % molSi/molN
% Biovars.RCN       = BioIn.RCN;              % molC/molN               
% Biovars.RCFe      = BioIn.RCFe;             % molC/molFe
Biovars.m0exp = zeros(nbsv,1);
Biovars.m0exp(1:Ewein.nlive)    = BioIn.m0exp;            % no unit

Biovars.Iopt                    = zeros(nbsv,1);
Biovars.Iopt(nemidx(1:11))      = Np.Iopt / 0.001433;     % W/m^2   
Biovars.alpha                   = zeros(nbsv,1);
Biovars.alpha(nemidx(1:11))     = Np.alpha * 0.001433;    % (W/m^2)^-1
Biovars.Vmax                    = zeros(nbsv,1);
Biovars.Vmax(nemidx(1:11))      = Np.Vmax;                % s^-1
Biovars.Kno3                    = zeros(nbsv,1);
Biovars.Kno3(nemidx(1:11))      = Np.Kno3 * 1000;         % molN/m^3
Biovars.pusai                   = zeros(nbsv,1);
Biovars.pusai(nemidx(1:11))     = Np.pusai / 1000;        % m^3/molN
Biovars.Knh4                    = zeros(nbsv,1);
Biovars.Knh4(nemidx(1:11))      = Np.Knh4 * 1000;         % molN/m^3 
Biovars.Ksi                     = zeros(nbsv,1);
Biovars.Ksi(nemidx(1:11))       = Np.Ksi * 1000;          % molSi/m^3
Biovars.Kgpp                    = zeros(nbsv,1);
Biovars.Kgpp(nemidx(1:11))      = Np.Kgpp;                % degC^-1
Biovars.gamma                   = zeros(nbsv,1);
Biovars.gamma(nemidx(1:11))     = Np.gamma;               % no unit
Biovars.res0                    = zeros(nbsv,1);
Biovars.res0(nemidx(1:11))      = Np.res0;                % s^-1
Biovars.Kres                    = zeros(nbsv,1);
Biovars.Kres(nemidx(1:11))      = Np.Kres;                % degC^-1
Biovars.Kdec                    = zeros(nbsv,nbsv);
Biovars.Kdec(nemidx(1:11),nemidx(1:11)) = Np.Kdec;        % degC^-1
Biovars.vdec                    = zeros(nbsv,nbsv);
Biovars.vdec(nemidx(1:11),nemidx(1:11)) = Np.vdec;        % s^-1

% Biovars.Kfe                     = zeros(nbsv,1);
% Biovars.Kfe(nemidx(1:2))        = BioIn.Kfe;              % molFe/m^3
% Biovars.kfe2n                   = zeros(nbsv,1);
% Biovars.kfe2n(nemidx(1:2))      = BioIn.kfe2n;            % molFe/molN
% Biovars.fe2nmax                 = zeros(nbsv,1);
% Biovars.fe2nmax(nemidx(1:2))    = BioIn.fe2nmax;          % molFe/molN
% Biovars.fe2nupfac               = BioIn.fe2nupfac;        % molFe/molN
% Biovars.ligbkg                  = BioIn.ligbkg;           % mol/m^3
% Biovars.alphascav               = BioIn.alphascav;        % s^-1
% Biovars.remineff                = BioIn.remineff;         % no unit
% Biovars.kiscav                  = BioIn.kiscav;           % W/m^2
% Biovars.kliglo                  = BioIn.kliglo;           % m^3/mol
% Biovars.klighi                  = BioIn.klighi;           % m^3/mol

% Iron parameters: Most scalars, a few assigned to phytoplankton groups
% only

Biovars.Kfe                     = zeros(nbsv,1);
Biovars.Kfe(nemidx(1:2))        = BioIn.Kfe;        % molFe/m^3
Biovars.kfe2n                   = zeros(nbsv,1);
Biovars.kfe2n(nemidx(1:2))      = BioIn.kfe2n;      % molFe/molN
Biovars.fe2nmax                 = zeros(nbsv,1);
Biovars.fe2nmax(nemidx(1:2))    = BioIn.fe2nmax;    % molFe/molN
Biovars.fe2nupfac               = BioIn.fe2nupfac;  % molFe/molN
Biovars.ligbkg                  = BioIn.ligbkg;     % mol/m^3
Biovars.alphascav               = BioIn.alphascav;  % s^-1
Biovars.remineff                = BioIn.remineff;   % no unit
Biovars.kliglo                  = BioIn.kliglo;     % mol lig^-1 kg^-1
Biovars.klighi                  = BioIn.klighi;     % mol lig^-1 kg^-1
Biovars.iofescav                = BioIn.iofescav;   % no unit

% Biovars.grmax     = Np.grmax;               % s^-1
% Biovars.lambda    = Np.lambda/1000;         % m^3/molN
% Biovars.thresh    = Np.thresh * 1000;       % molN/m^3
% Biovars.Kgra      = Np.Kgra;                % degC^-1
% Biovars.mor0      = Np.mor0/1000;           % (molN/m^3)^-1 s^-1
% Biovars.Kmor      = Np.Kmor;                % degC^-1
% Biovars.vdec      = Np.vdec;                % s^-1
% Biovars.Kdec      = Np.Kdec;                % degC^-1


% Nekton: the paramters for nekton come from the Ecopath mass balance

Biovars.b0 = zeros(nbsv,1);
Biovars.b0(1:Ewein.ngroup) = Ep.b;
Biovars.q0 = zeros(nbsv);
Biovars.q0(1:Ewein.nlive,1:Ewein.nlive) = Ep.q0(1:Ewein.nlive,1:Ewein.nlive);
Biovars.ge = zeros(nbsv,1);
Biovars.ge(1:Ewein.nlive) = Ep.ge(1:Ewein.nlive);
Biovars.gs = zeros(nbsv,1);
Biovars.gs(1:Ewein.nlive) = Ewein.gs(1:Ewein.nlive);


% For x, d, and theta, I accept data in one of two ways:
% 1) vector of P-value log-tranformed-anomaly-from-base values, identical
% to those used in the functional response file input for aydin-ecosim
% 2) matrix of values for each group pair.  These values are NOT anomalies
% but actual values to be used in the functional response equations.

if isvector(Ewein.x)
    [xj,xi] = meshgrid(Ewein.x);
    x = 1 + exp(xi + xj);
    x(Ewein.dc == 0) = 0;
else
    x = Ewein.x;
end
    
if isvector(Ewein.d)
    [dj,di] = meshgrid(Ewein.d);
    d = 1 + exp(di + dj);
    d(Ewein.dc == 0) = 0;
else
    d = Ewein.d;
end

if isvector(Ewein.theta)
    [thj,thi] = meshgrid(Ewein.theta);
    theta = exp(0.05*(thi + thj)); % TODO double-check theta, since different in spreadsheet and ppt
    theta(Ewein.dc == 0) = 0;
else
    theta = Ewein.theta;
end
    
n = nbsv - Ewein.ngroup;
Biovars.x = padarray(x, [n n], 'post');
Biovars.d = padarray(d, [n n], 'post');
Biovars.theta = padarray(theta, [n n], 'post');

% Zooplankton: grazing same as predation, but in volumetric terms

dz = mld; % TODO: This is arbitrary, and possibly very important... run some tests

Biovars.b0v = Biovars.b0./dz; % mol N m^-3
Biovars.q0v = Biovars.q0./dz; % mol N m^-3 s^-1

% Production (for debugging with ecosim production only)

Biovars.p0 = zeros(nbsv,1);
Biovars.p0(1:Ewein.ngroup) = Ep.pb .* Ep.b; % mol N m^-2 s^-1
Biovars.p0v = Biovars.p0./dz;                     % mol N m^-3 s^-1 

% Temperature factors

Biovars.Kgra = zeros(nbsv,1);
Biovars.Kgra(iszoo & ~isextrazoo) = Np.Kgra(3:5);
Biovars.Kgra(isextrazoo) = BioIn.kgra;

tempfac0 = exp(Biovars.Kgra .* BioIn.temp);

Biovars.q0vat0 = bsxfun(@rdivide, Biovars.q0v, tempfac0');

% Mortality (using ecopath-based mortality)

Biovars.m0 = zeros(nbsv,1);
Biovars.m0(1:Ewein.nlive) = Ep.otherMortRate(1:Ewein.nlive);  % s^-1
if any(Biovars.m0 < 0)
    warning('WCE:m0neg', 'Some m0 < 0, setting to 0');
    Biovars.m0 = max(Biovars.m0, 0);
end

% Mortality: Allows for any exponential function aB.^y

m0 = zeros(nbsv,1);
m0(1:Ewein.nlive) = Ep.otherMortRate(1:Ewein.nlive);  % s^-1
m0b = m0 .* bio(1,:)'; % mass-balanced flux, molN/m^3/s

Biovars.m0coef = m0b./(bio(1,:)'.^Biovars.m0exp);
Biovars.m0coef = max(Biovars.m0coef, 0); % get rid of NaNs (and I guess negatives, if they ever appeared, though they shouldn't)

if Biovars.diapause
    bzl = sum(bio(1,[Biovars.idx.zl1 Biovars.idx.zl2]));
    m0bzl = Ep.otherMortRate(Biovars.idx.zl) .* bzl;
    mexp = Biovars.m0exp(Biovars.idx.zl);
    Biovars.m0coef(Biovars.idx.zl) = 0;
    Biovars.m0coef([Biovars.idx.zl1 Biovars.idx.zl2]) = m0bzl./(bzl.^mexp);
end

% Biovars.m0max = zeros(nbsv,1);
% Biovars.m0max(1:Ewein.nlive) = Ep.otherMortRate(1:Ewein.nlive);  % s^-1
% Biovars.m0max = max(Biovars.m0max,0);
% 
% Biovars.m0coef = Biovars.m0max./(bio(1,:)'.^Biovars.m0exp);
% Biovars.m0coef = max(Biovars.m0coef, 0); % get rid of NaNs (and I guess negatives, if they ever appeared, though they shouldn't)

% Biovars.dm0db = Biovars.m0max./bio(1,:)'; % s^-1 (mol N m ^-3)^-1
% Biovars.dm0db = max(Biovars.dm0db,0); % get rid of NaNs (and I guess negatives, if they ever appeared, though they shouldn't)

% PON and Opal sink

Biovars.settle = zeros(nz,nbsv);
Biovars.settle(:,nemidx(1:11)) = repmat(-Np.settle', nz, 1);    % m/s
Biovars.settle(:,Biovars.idx.pofe) = Biovars.settle(:,Biovars.idx.pon); % Same as PON

% A few odds and ends

Biovars.ecosimppflag = BioIn.ecosimpp;

% Diapause setup


if Biovars.diapause
    
    epvars = {'gs', 'ge', 'm0exp'};
    for ii = 1:length(epvars)
        Biovars.(epvars{ii})([Biovars.idx.zl1 Biovars.idx.zl2]) = Biovars.(epvars{ii})(Biovars.idx.zl);
    end
    
    % Set when they swim up and down
    
    dv = datevec(datenum(Grd.start_date) + Grd.time/86400);
    day = datenum(dv) - datenum([dv(:,1) ones(Grd.nt,2)]);
    
    Biovars.zlswim = zeros(Grd.nt,1);
    Biovars.zlswim(day >= BioIn.dday(1) & day <= BioIn.dday(2)) = 1; % swim up
    Biovars.zlswim(day >= BioIn.dday(3) | day <  BioIn.dday(1)) = -1; % stay down
    
    % Set when to split and recombine the two groups
    
    Biovars.zlsplit   = [false; Biovars.zlswim(2:end) == -1 & Biovars.zlswim(1:end-1) ~= -1];
    d1 = find(Biovars.zlsplit, 1);
    Biovars.zlswim(1:d1-1) = 0; % keep as one until first down day
    Biovars.zlcombine = [false; Biovars.zlswim(2:end) == 0 & Biovars.zlswim(1:end-1) == 1];
    
    % Some extras
    
    Biovars.t = Grd.time;
    Biovars.dfrac = BioIn.dfrac;
    Biovars.eatdcop = BioIn.eatdcop;
    
end

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt)

Param       = Biovars;
Param.irr   = meanqi; 
Param.temp  = temperature;
Param.z     = z;
Param.dz    = dz;
    
if any(isnan(oldbio(:)))
    warning('WCE: NaN in biology');
end

% Integrate biology over this time step

if Biovars.diapause
    it = find(t == Biovars.t);
    zltot = sum(oldbio(:,[Biovars.idx.zl1 Biovars.idx.zl2]), 2);
    if Biovars.zlsplit(it)
        oldbio(:,Biovars.idx.zl2) = Biovars.dfrac .* zltot;
        oldbio(:,Biovars.idx.zl1) = zltot - oldbio(:,Biovars.idx.zl2);
    elseif Biovars.zlcombine(it)
        oldbio(:,Biovars.idx.zl1) = zltot;
        oldbio(:,Biovars.idx.zl2) = 0;
    end
    oldbio(:,Biovars.idx.zl) = 0;
end

[newbio, db, Flx, Diag, badthings] = integratebio(@wceode, t, dt, oldbio, Param, Biovars.odesolver{:});

if Biovars.diapause
    newbio(:,Biovars.idx.zl) = newbio(:,Biovars.idx.zl1) + newbio(:,Biovars.idx.zl2);
end

% Check and correct for silica issue

isneg = newbio(:, Biovars.idx.plsi) < 0;
Diag.extrasi = zeros(size(isneg));
Diag.extrasi(isneg) = -newbio(isneg,Biovars.idx.plsi);
newbio(isneg, Biovars.idx.plsi) = 0;
badthings(isneg, Biovars.idx.plsi) = false;

% Check for problems

if any(badthings(:))
    [ridx,cidx] = find(badthings);
    nb = length(ridx);
    errstr = cell(nb,1);
    for ii = 1:nb
        badbio = newbio(ridx(ii), cidx(ii));
        if isnan(badbio)
            errstr{ii} = sprintf('NaN: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), t);
        elseif isinf(badbio)
            errstr{ii} = sprintf('Inf: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), t);
        elseif badbio < 0
            errstr{ii} = sprintf('Neg: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), t); 
        end
    end
    errstr = sprintf('  %s\n', errstr{:});
    errstr = sprintf('Biology out of range:\n%s', errstr);
    error('WCE:biologyOutOfRange', errstr);
end

% Diagnostics: Intermediate fluxes

nd = size(Biovars.wceflux,1);

fluxes = zeros(length(z),nd);
for id = 1:nd
    fluxes(:,id) = Flx(1).(Biovars.wceflux{id,1})(Biovars.wceflux{id,2},Biovars.wceflux{id,3},:);
end

% Diagnistics: Other

ndiag = 15;
if isempty(Diag) % If use any solver other than euler
    diag = zeros(length(z), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(2, diag{:});
end

diag = [db fluxes diag];

% Conservation check (debugging, very ESA-specific)

% tot = sum(newbio .* Param.dz, 1); 
% ntot = sum(tot(1:27));          % mol m^-2
% stot = sum(tot([28:29 31]));    % mol m^-2
% ftot = sum(tot([30 32:33]));    % umol m^-2
% fprintf(Biovars.cfid, '%f %f %f %f\n', t, ntot, stot, ftot);
% if abs(ntot - Biovars.ntot0) > 0.01
%     fclose(Biovars.cfid);
%     error('WCE:notConserved', 'N not conserved');
% end
                      
%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = Biovars.settle;

if Biovars.diapause
    
    it = find(Biovars.t == t);
    
    if Biovars.zlswim(it) == 1
        wsink(z < -10, Biovars.idx.zl2) = 80./86400; % swim up
    elseif Biovars.zlswim(it) == -1
        wsink(z > -400, Biovars.idx.zl2) = -80./86400; % stay down
        wsink(z < -450, Biovars.idx.zl2) = 80./86400;
    end
    
end

