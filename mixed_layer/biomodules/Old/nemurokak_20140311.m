function varargout = nemurokak(action, varargin)
%NEMUROKAK NEMURO, with iron and grazing modifications
%
% See biomodule.m for full syntax details.
%
% This module runs a lower-trophic level biogeochemical model based on the
% NEMURO model.  It includes modifications for explicit iron limitation, as
% well as options to switch between single-resource and multi-resource
% grazing functional responses.
%
% User-specified input variables (passed to mixed_layer as parameter/value
% pairs):
%
%   NemParam:   1 x 1 structure of nemuro input variables (see
%               nemuroinputparser.m)
%
%   bnem0:      nz x 13 array of initial values for each state variables.
%               Column 1 holds depths (m, negative down) and columns 2-13
%               correspond to the 12 state variables (mol/m^3):
%               1: PS   4: ZL   7: NH4  10: SiOH4
%               2: PL   5: ZP   8: PON  11: Opal
%               3: ZS   6: NO3  9: DON  12: Fe
%
%   ivlev:      String, specifying which grazing scheme to use ['orig']
%               'orig':     single-resource for all except ZP, which has
%                           "gourmet" function
%               'single':   single-resource for all  
%               'multi':    multi-resource for all
%
%   grmax:      (ivlev = 'multi' only) 11 x 1 array, maximum razing rates
%               for each of the live state variables (with placeholders for
%               the remaining originals) (/s)
%
%   thresh:     (ivlev = 'multi' only) 11 x 1 array, theshold feeding
%               values for each of the live state variables (mol/m^3)
%
%   p:          (ivlev = 'multi' only) 11 x 11 array of prey preference
%               values, with values ranging from 0-1.  1 indicates the most
%               preferred prey.
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
%   diapause:   logical scalar, if true, implement seasonal migration for
%               the ZL group
%
%   New iron
%   remineff:   remineralization efficiency for sinking iron detritus
%               relative to organic material [0.5] 
%
%   kliglo:     ligand binding coefficient in low light conditions (mol
%               lig^-1 kg^-1) [1.0e12]
%
%   klighi:     ligand binding coefficient in high light conditions (mol
%               lig^-1 kg^-1) [1.0e8]
%
%   iofescav:   light threshold for transition between weak and strong iron
%               binding (W m^-2) [10]
%
%   alphascav:  iron scavenging coefficient onto sinking detritus (s^-1)
%               [15.0/(86400*365)]
%
%   felig2don:  coefficient for scaling ligand concentration with organic
%               nitrogen (mol lig (mol N^-1)) [0]
%
%   ligbkg:     background ligand concentration (mol lig kg^-1) [1e-9]


% Copyright 2011 Kelly Kearney

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


p.addParamValue('NemParam',  [],                        @isstruct);
p.addParamValue('bnem0',     []);
p.addParamValue('ivlev',     'orig',                    @(x) ischar(x) && ismember(x, {'single', 'multi', 'orig', 'mishmash'}));
% p.addParamValue('RCN',       6.625,               @(x) isnumeric(x) && iscalar(x));
% p.addParamValue('RCFe',      33e3,                @(x) isnumeric(x) && iscalar(x));
% p.addParamValue('Kfe',       [0.033e-6 0.188e-6], @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('odesolver', {'euler'});
p.addParamValue('reroute',   cell(0,5),                 @(x) isempty(x) || (iscell(x) && size(x,2)==5));
p.addParamValue('grmax',     [],                        @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 11));
p.addParamValue('thresh',    [],                        @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 11));
p.addParamValue('p',         [],                        @(x) isempty(x) || (isnumeric(x) && isequal(size(x), [11 11])));
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
p.addParamValue('diapause',  false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('dday',      [120 134 243],             @(x) isnumeric(x) && isvector(x) && length(x)==3);
p.addParamValue('dfrac',     0.9,                       @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);


p.parse(In);
BioIn = p.Results;

nodefault = {'NemParam','bnem0'};
tf = ismember(p.UsingDefaults, nodefault);
missing = sprintf('%s, ', p.UsingDefaults{tf});
if any(tf)
    error('Missing input for nemurokak module: %s', missing(1:end-2));
end

%------------------------------
% Initial biomass
%------------------------------

S = warning('off', 'MATLAB:interp1:NaNinY');

% Check that bnem0 has correct # of columns

nbsv = 18; % Number of tracked biological state variables

nz = length(Grd.z);
nb0 = size(BioIn.bnem0,2);
if nb0 ~= 13
    error('Initial biomass (bnem0) should include depth plus 12 state variables');
end

bio = zeros(nz, nbsv);

% Most state variables' initial profiles set by input...

if all(BioIn.bnem0 >= 0)
    BioIn.bnem0(:,1) = -BioIn.bnem0(:,1); % Tired of trying to remember +up or down
end

bio(:,1:12) = interp1(BioIn.bnem0(:,1), BioIn.bnem0(:,2:end), Grd.z); % mol N/m^3

% ... except PL silica, which is proportional to N

bio(:,13) = bio(:,2) .* BioIn.NemParam.RSiN;

% ... and PS and PL iron, which starts at 0.  TODO: Or should I assume
% something else?  Half-sat Fe:N maybe?

% bio(:,14:15) = bsxfun(@times, bio(:,1:2), BioIn.fe2nmax);
bio(:,14:15) = 0;

% ... and POFe starts at 0,

bio(:,16) = 0;

% If diapause option is on, split ZL into two groups.  Start with all
% biomass in the non-migrating group.

if BioIn.diapause
    zltot = bio(:,4);
    bio(:,17) = zltot;
    bio(:,18) = 0;
    bio(:,4)  = 0;
end

warning(S);

%----------------------
% Mixing and forcing
%----------------------

% All variables mixed

ismixed = true(1, nbsv);

% No bottom forcing

bottomval = nan(1, nbsv);

%----------------------
% Variable names
%----------------------

% State variables

names = {
    'PS',       'Small Phytoplankton',          'molN/m^3'
    'PL',       'Large Phytoplankton',          'molN/m^3'
    'ZS',       'Small Zooplankton',            'molN/m^3'
    'ZL',       'Large Zooplankton',            'molN/m^3'
    'ZP',       'Pradatory Zooplankton',        'molN/m^3'
    'NO3',      'Nitrate',                      'molN/m^3'
    'NH4',      'Ammmonium',                    'molN/m^3'
    'PON',      'Particulate Organic Nitrogen', 'molN/m^3'
    'DON',      'dissolved Organic Nitrogen',   'molN/m^3'
    'SiOH4',    'Silicate',                     'molSi/m^3'
    'Opal',     'Particulate Opal',             'molSi/m^3'
    'Fe'        'Dissolved Iron'                'molFe/m^3'
    'PLsi'      'Large Phytoplankton Silica',   'molSi/m^3'
    'PSfe'      'Small Phytoplankton Iron'      'molFe/m^3'
    'PLfe'      'Large Phytoplankton Iron'      'molFe/m^3'
    'POFe'      'Particulate iron'              'molFe/m^3'
    'ZL1'       'Large Zooplankton (no diapause)' 'molN/m3'
    'ZL2'       'Large Zooplankton (diapause)'  'molN/m^3'};

% Indices

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
Biovars.idx.psfe  = find(strcmp(names(:,1), 'PSfe'));
Biovars.idx.plfe  = find(strcmp(names(:,1), 'PLfe'));
Biovars.idx.mys   = nbsv + 1; % "Mystery" box, for source/sink terms coming from or going to nowhere
Biovars.idx.zs    = find(strcmp(names(:,1), 'ZS'));
Biovars.idx.zl    = find(strcmp(names(:,1), 'ZL'));
Biovars.idx.zp    = find(strcmp(names(:,1), 'ZP'));
Biovars.idx.pofe  = find(strcmp(names(:,1), 'POFe'));
Biovars.idx.zl1   = find(strcmp(names(:,1), 'ZL1'));
Biovars.idx.zl2   = find(strcmp(names(:,1), 'ZL2'));


% Diagnostics

diagnames = {...
    'PSpsmax',      'Temp-mediated photsynthesis (small)'    '/s'
    'PLpsmax',      'Temp-mediated photsynthesis (large)'    '/s'
    'PSlightlim',   'Light limitation (small)',         'no units'
    'PLlightlim',   'Light limitation (large)',         'no units'
    'PSno3lim',     'Nitrate limitation (small)',       'no units'
    'PLno3lim',     'Nitrate limitation (large)',       'no units'
    'PSnh4lim',     'Ammonium limitation (small)',      'no units'
    'PLnh4lim',     'Ammonium limitation (large)',      'no units'
    'PLsilim',      'Silica limitation (large)'         'no units'
    'PSfelim'       'Iron limitation (small)'           'no units'
    'PLfelim'       'Iron limitation (large)'           'no units'
    'PStemplim'     'Temperature factor (small)',       'no unit'
    'PLtemplim'     'Temperature factor (large)',       'no unit'
    'I',            'Irradiance'                        'W m^-2'
    'kappa',        'Attenuation coefficient'           'm^-1'
    'kp'            'Attenuation self-shading only',    'm^-1'
    'PSfe2n'        'Fe:N ratio (small)'                'molFe/molN'
    'PLfe2n'        'Fe:N ratio (large)'                'molFe/molN'
    'PSfedef'       'Fe deficiency (small)'             'molFe/molN'
    'PLfedef'       'Fe deficiency (large)'             'molFe/molN'};

% Intermediate fluxes

reroute = BioIn.reroute;
[tf, loc] = ismember(reroute(:,2:4), names);
if ~all(tf)
    error('Unrecognized critter in reroute table');
end
reroute(:,2:4) = num2cell(loc);

fluxlist = listfluxes('nemurokak', Biovars.idx);

fluxlist = [fluxlist; reroute(:,[1 2 4])];
nd = size(fluxlist,1);

fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', fluxlist{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', fluxlist{id,:});
    fluxnames{id,3} = 'mol/m^3/s'; 
end

Biovars.nemflux = fluxlist;

% Combine intermediate fluxes and extra diagnostics

diagnames = [diagnames; fluxnames];

%------------------------------
% Variables needed for ODE
%------------------------------

% ODE solver

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

% Rerouting of fluxes

Biovars.reroute = reroute;

% Extended NEMURO parameters

Np = nemuroflexinput(BioIn.NemParam);

Biovars.alpha1    = Np.alpha1;              % m^-1
Biovars.alpha2    = Np.alpha2/1000;         % (m^3/molN)/m
Biovars.usesteele = Np.usesteele;           % logical
Biovars.Iopt      = Np.Iopt / 0.001433;     % W/m^2        
Biovars.alpha     = Np.alpha * 0.001433;    % (W/m^2)^-1
Biovars.Vmax      = Np.Vmax;                % s^-1
Biovars.Kno3      = Np.Kno3 * 1000;         % molN/m^3
Biovars.pusai     = Np.pusai / 1000;        % m^3/molN
Biovars.Knh4      = Np.Knh4 * 1000;         % molN/m^3 
Biovars.Ksi       = Np.Ksi * 1000;          % molSi/m^3
Biovars.Kgpp      = Np.Kgpp;                % degC^-1
Biovars.RSiN      = Np.RSiN;                % molSi/molN
Biovars.gamma     = Np.gamma;               % no unit
Biovars.res0      = Np.res0;                % s^-1
Biovars.Kres      = Np.Kres;                % degC^-1
Biovars.grmax     = Np.grmax;               % s^-1
Biovars.lambda    = Np.lambda/1000;         % m^3/molN
Biovars.thresh    = Np.thresh * 1000;       % molN/m^3
Biovars.Kgra      = Np.Kgra;                % degC^-1
Biovars.mor0      = Np.mor0/1000;           % (molN/m^3)^-1 s^-1
Biovars.Kmor      = Np.Kmor;                % degC^-1
Biovars.vdec      = Np.vdec;                % s^-1
Biovars.Kdec      = Np.Kdec;                % degC^-1

switch BioIn.ivlev
    case 'orig'

        % Need these ones for originl Ivlev only

        Biovars.psiPL = BioIn.NemParam.PusaiPL/1000; % m^3/molN
        Biovars.psiZS = BioIn.NemParam.PusaiZS/1000; % m^3/molN
        
    case 'multi'
        
        % User-input matrices overshadow NP set
        
        if isempty(BioIn.grmax) || isempty(BioIn.thresh) || isempty(BioIn.p)
            error('Need to provide grmax, thresh, and p for multi-resource');
        end
        
        Biovars.grmax  = BioIn.grmax;
        Biovars.thresh = BioIn.thresh;
        Biovars.p      = BioIn.p;
        
    case 'mishmash'
        
        Biovars.psiPL = BioIn.NemParam.PusaiPL/1000; % m^3/molN
        Biovars.psiZS = BioIn.NemParam.PusaiZS/1000; % m^3/molN
        if isempty(BioIn.grmax) || isempty(BioIn.thresh) || isempty(BioIn.p)
            error('Need to provide grmax, thresh, and p for multi-resource');
        end
        
        Biovars.grmax  = BioIn.grmax;
        Biovars.thresh = BioIn.thresh;
        Biovars.p      = BioIn.p;
        
end
     
Biovars.ivlev          = BioIn.ivlev;    
% Biovars.RCN       = BioIn.RCN;            % molC/molN               
% Biovars.RCFe      = BioIn.RCFe;           % molC/molFe
Biovars.gs              = 1 - Np.alphaeg;   % no unit
Biovars.ge              = Np.beta;          % no unit

% Iron-related

Biovars.Kfe             = zeros(nbsv,1);
Biovars.Kfe(1:2)        = BioIn.Kfe;        % molFe/m^3
Biovars.kfe2n           = zeros(nbsv,1);
Biovars.kfe2n(1:2)      = BioIn.kfe2n;      % molFe/molN
Biovars.fe2nmax         = zeros(nbsv,1);
Biovars.fe2nmax(1:2)    = BioIn.fe2nmax;    % molFe/molN
Biovars.fe2nupfac       = BioIn.fe2nupfac;  % molFe/molN
Biovars.ligbkg          = BioIn.ligbkg;     % mol/m^3
Biovars.alphascav       = BioIn.alphascav;  % s^-1
Biovars.remineff        = BioIn.remineff;   % no unit
Biovars.kliglo          = BioIn.kliglo;     % mol lig^-1 kg^-1
Biovars.klighi          = BioIn.klighi;     % mol lig^-1 kg^-1
Biovars.iofescav        = BioIn.iofescav;   % no unit

% Diapause-related

Biovars.diapause        = BioIn.diapause;
if Biovars.diapause
    dv = datevec(datenum(Grd.start_date) + Grd.time/86400);
    day = datenum(dv) - datenum([dv(:,1) ones(Grd.nt,2)]);
    
    % Set when they swim up and down
    
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
end

Biovars.ge(Biovars.idx.zl1)   = Biovars.ge(Biovars.idx.zl);
Biovars.gs(Biovars.idx.zl1)   = Biovars.gs(Biovars.idx.zl);
Biovars.mor0(Biovars.idx.zl1) = Biovars.mor0(Biovars.idx.zl);
Biovars.Kmor(Biovars.idx.zl1) = Biovars.Kmor(Biovars.idx.zl);
Biovars.ge(Biovars.idx.zl2)   = Biovars.ge(Biovars.idx.zl);
Biovars.gs(Biovars.idx.zl2)   = Biovars.gs(Biovars.idx.zl);
Biovars.mor0(Biovars.idx.zl2) = Biovars.mor0(Biovars.idx.zl);
Biovars.Kmor(Biovars.idx.zl2) = Biovars.Kmor(Biovars.idx.zl);

% PON, Opal, and POFe sink

Biovars.sink = zeros(1,nbsv);
Biovars.sink(1,Biovars.idx.pon)  = -Np.settle(Biovars.idx.pon); 
Biovars.sink(1,Biovars.idx.opal) = -Np.settle(Biovars.idx.opal);
Biovars.sink(1,Biovars.idx.pofe) = -Np.settle(Biovars.idx.pon); % Same as PON

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
Param       = Biovars;
Param.irr   = meanqi; 
Param.temp  = temperature;
Param.z     = z;
Param.dz    = dz;
             
% For diapause splitting/combining, need to do this outside the ODE

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

[newbio, db, Flx, Diag, badthings] = integratebio(@nemurokakode, t, dt, oldbio, Param, Biovars.odesolver{:});

if Biovars.diapause
    newbio(:,Biovars.idx.zl) = newbio(:,Biovars.idx.zl1) + newbio(:,Biovars.idx.zl2);
end

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
    error('NEMURO:biologyOutOfRange', errstr);

end

% If diapause, add the ZL groups

% Diagnostics

ndiag = 14;

if isempty(Diag) % If use any solver other than euler
    diag = zeros(length(z), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(2, diag{:});
end

% nemflux = findnemflux(nemuroflexinput(Biovars), 'list');
% nemflux{1} = [nemflux{1} Biovars.reroute(:,1)'];            % add rerouted
% nemflux{2} = [nemflux{2} cell2mat(Biovars.reroute(:,2)')];  % add rerouted
% nemflux{3} = [nemflux{3} cell2mat(Biovars.reroute(:,4)')];  % add rerouted

nemflux = Biovars.nemflux;
nfx = size(nemflux,1);
nz = size(newbio,1);
fluxes = zeros(nz, nfx);
for ifx = 1:nfx
    fluxes(:,ifx) = Flx(1).(nemflux{ifx,1})(nemflux{ifx,2}, nemflux{ifx,3},:);
end

diag = [diag fluxes];

%**************************************************************************


function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = ones(size(z)) * Biovars.sink;

if Biovars.diapause
    it = find(Biovars.t == t);
    
    if Biovars.zlswim(it) == 1
        wsink(z < -10, Biovars.idx.zl2) = 80./86400; % swim up
    elseif Biovars.zlswim(it) == -1
        wsink(z > -400, Biovars.idx.zl2) = -80./86400; % stay down
        wsink(z < -450, Biovars.idx.zl2) = 80./86400;
    end
    
end




