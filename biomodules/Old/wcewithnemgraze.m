function varargout = wcewithnemgraze(action, varargin)
%WCE Water column ecosystem biological module
%
% See biomodule.m for full syntax details.
%
% This module simulates a mixed planktonic-nektonic ecosystem, based on
% a combination of models derived mainly from NEMURO and Kerim Aydin's
% version of Ecosim.
%
% A note on units: The biomass of all critters is saved to file in mol
% N/m^3 (or umol Fe/m^3).  For nektonic critters, all biomass is placed in
% the surface cell, and actually represents the total over the entire water
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
%   reroute:    n x 5 cell array.  This allows you to reroute fluxes from
%               the original path defined in the nemuro model.  Each row
%               desribes as change, with column as follows:
%               col 1:  name of flux (gpp, graze, resp, excrete, exc,
%                       egest, mort, dec, ferelax)
%               col 2:  index of original source group
%               col 3:  index of original sink group
%               col 4:  index of new sink group
%               col 5:  fraction of flux to reroute
%
%   Fe:         1 x 1 structure of iron-related parameters. The structure
%               holds the following fields:
%   
%               tfe:    timescale for iron relaxation (s)
%
%               a:      empirical Fe:C power
%
%               b:      empirical Fe:C coefficient (molC m^-3)^-1, ...units
%                       are a bit fuzzy, based on the power law, but this
%                       is how Fiechter defines them)
%
%               kfe:    Fe:C half-saturation constant (umolFe/molC)
%
%               frem:   fraction of iron remineralized 
%
%               cpern:  C:N ratio (molC/molN)           
%
%               fe0:    n x 2 array defining initial profile of iron.
%                       Column 1 holds depth (m), column 2 iron
%                       concentration (umol/m^3)
%
%   odesolver:  cell array of strings, indicating which solvers to use.  If
%               the first one fails to integrate a timestep (i.e. causes
%               something to become negative, NaN, or Inf), the next one is
%               tried.  See integratebio.m for choices. [{'euler'}]
%
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
p.addParamValue('Fe', [], @isstruct);

p.parse(In);
BioIn = p.Results;

nodefault = {'Ewein', 'NemParam', 'types', 'bnem0', 'Fe'};
tf = ismember(p.UsingDefaults, nodefault);
missing = sprintf('%s, ', p.UsingDefaults{tf});
if any(tf)
    error('Missing input for wce module: %s', missing(1:end-2));
end

%------------------------------
% Set up biology
%------------------------------

[bio, names, ismixed, bottomval, Biovars] = wcesetup(...
    Grd, BioIn.Ewein, BioIn.NemParam, BioIn.types, BioIn.bnem0, ...
    BioIn.nomix, BioIn.mld, BioIn.temp, BioIn.kgra, BioIn.Fe);

% Add flag for ecosim primary production (mostly for debugging)

Biovars.ecosimppflag = BioIn.ecosimpp;

%------------------------------
% Ode solver
%------------------------------

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

%------------------------------
% Diagnostics
%------------------------------

% db/dt for each group

dbnames = names;
for idb = 1:size(dbnames,1)
    dbnames{idb,1} = ['d' dbnames{idb,1}];
    dbnames{idb,2} = ['dB/dt: ' dbnames{idb,2}];
    dbnames{idb,3} = 'molN m^-3 s^-1';
end

% Fluxes between groups, by process type

diaglist = fluxindices('wce', Biovars.links, Biovars.idx, Biovars.m0, Biovars.vdec);
if Biovars.ecosimppflag
    ismissing = strcmp(diaglist(:,1), 'gpp') | ...
                strcmp(diaglist(:,1), 'exc') | ...
                strcmp(diaglist(:,1), 'resp');
else
    ismissing = strcmp(diaglist(:,1), 'npp');
end
diaglist = diaglist(~ismissing,:);

% Add user-rerouted fluxes

if isempty(BioIn.reroute)
    Biovars.reroute = cell(0,5);
else
    Biovars.reroute = BioIn.reroute;
    [tf, idx] = ismember(Biovars.reroute(:,2:4), names(:,2));
    Biovars.reroute(:,2:4) = num2cell(idx);
end

diaglist = [diaglist; Biovars.reroute(:,[1 2 4])];
nd = size(diaglist,1);

diagnames = cell(nd,3);
for id = 1:nd
    diagnames{id,1} = sprintf('%s_%02d_%02d', diaglist{id,:});
    diagnames{id,2} = sprintf('%s: %02d to %02d', diaglist{id,:});
    diagnames{id,3} = 'mol N m^-3 s^-1';
end

Biovars.diaglist = diaglist;

% Some intermediates related to primary production

diagnames2 = {
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
    'extrasi'       'Extra silica due to negatives'     'mol Si m^-3'};

diagnames = [dbnames; diagnames; diagnames2];

% Additional stuff for wce variants

Biovars.Np = BioIn.NemParam;

Biovars.idx.zs = find(strcmp(names(:,1), 'ZS'));
Biovars.idx.zl = find(strcmp(names(:,1), 'ZL'));
Biovars.idx.zp = find(strcmp(names(:,1), 'ZP'));

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);

Param       = Biovars;
Param.irr   = meanqi; 
Param.temp  = temperature;
Param.z     = z;
Param.dz    = dz;
    
if any(isnan(oldbio(:)))
    warning('WCE: NaN in biology');
end

% Integrate biology over this time step

[newbio, db, WceOut, Diag, badthings] = integratebio(@wcewithnemgrazeode, t, dt, oldbio, Param, Biovars.odesolver{:});

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

nd = size(Biovars.diaglist,1);

diag = zeros(length(z),nd);
for id = 1:nd
    diag(:,id) = WceOut(1).(Biovars.diaglist{id,1})(Biovars.diaglist{id,2},Biovars.diaglist{id,3},:);
end

% Diagnistics: Other

ndiag = 12;
if isempty(Diag) % If use any solver other than euler
    diag2 = zeros(length(z), ndiag);
else
    diag2 = struct2cell(Diag);
    diag2 = cat(2, diag2{:});
end


diag = [db diag diag2];

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

% MODIFY CODE HERE

wsink = Biovars.settle;