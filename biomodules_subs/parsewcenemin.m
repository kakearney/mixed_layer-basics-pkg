function BioIn = parsewcenemin(varargin)
%PARSEWCWNEMIN Parse input for wce or nemurokak mixed_layer modules
%
% Optional input variables (passed as parameter/value pairs)
%
% Applicable to both nemurokak and wce:
%
%   NemParam:   1 x 1 structure of NEMURO parameters (see
%               nemuroinputparser) 
%
%   bnem0:      nz x 13 array of initial values for each state variable.
%               Column 1 holds depths (m, negative down) and columns 2-14
%               correspond to the 13 state variables (mol/m^3):
%               1: PS   4: ZL   7: NH4  10: SiOH4  13: POFe
%               2: PL   5: ZP   8: PON  11: Opal
%               3: ZS   6: NO3  9: DON  12: Fe 
%
%   odesolver:  cell array of strings, indicating which solvers to use.  If
%               the first one fails to integrate a timestep (i.e. causes
%               something to become negative, NaN, or Inf), the next one is
%               tried.  See integratebio.m for choices. [{'euler'}]
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
%   remineff:   Fraction of particulate iron remineralization relative to
%               organic (N) material (no unit) [0.25]
%
%   kliglo:     Lower limit of ligand binding under low-light conditions
%               (m^3/mol) [3.0e8]
%
%   klighi:     Upper limit of ligand binding under high-light conditions
%               (m^3/mol) [1.0e5]  
%
%   iofescav:   Factor by which ligand-binding strength decreases with
%               increasing light levels 
%
%   ligbkg:     Ligand background concentration (mol/m^3) [1.0e-6]
%
%   kfe2n:      1 x 2 array, half-saturation constants for internal Fe:N
%               ratio for small and large phytoplankton, respectively
%               (molFe/molN) [[6.625e-05 0.0001325], i.e. [10 20]
%               umolFe/molC assuming Redfield]
%
%   Kfe:        1 x 2 array, half-saturation constants for iron uptake for
%               small and large phytoplankton, respectively (molFe/m^3)
%               [6e-7 3e-6]
%
%   fe2nmax:    1 x 2 array, maximum internal Fe:N ratio for small and
%               large phytoplankton, respectively (molFe/molN) 
%               [[0.00033125 0.0033125], i.e. [50 500] umolFe/molC assuming
%               Redfield]
%
%   fe2nupfac:  scalar, Fe:N uptake ratio (molFe/molN) [15e-6]
%
%   alphascav:  Iron scavenging coefficient (s^-1) [4.7565e-07, i.e. 15
%               yr^-1]
%
%   scavthresh: level of dissolved iron at which scavenging onto particles
%               begins to increase (mol Fe/m^3) [1e-6]
%
%   scavfac:    factor by which the scavenging rate increases once over the
%               above threshold ((molFe/m^3)^-1) [0]
%
%   diapause:   logical scalar, true to turn on diapause for ZL group
%               [false]
%
%   dday:       1 x 3 vector, day of year when diapausing copepods begin
%               migrating up, stop migrating up, and begin migrating down,
%               respectively [120 134 243]
%
%   dfrac:      scalar, fraction of copepod population that migrates [0.9]
%
%   preyvis:    n x m array idicating how much of the prey fields each
%               predator can see.  The first column holds depth values, and
%               must span the full water column.  The remaining columns
%               correspond to all state variables in the model, an indicate
%               fraction of the water column "seen" (0-1).  By default, all
%               groups see the full water column.  This is most relevant to
%               nektonic predators in wce, but I've left the option for
%               plankton as well.
%
% WCE only:
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
%   m0exp:      Exponent for mortality function, of form M0 = aB^(m0exp). A
%               value of 1 leads to linear mortality and 2 to quadratic
%               mortality. Can also be a nlive x 1 array of values to allow
%               different functions for each critter.[2] 
%
% NEMUROKAK Only:
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

% Copyright 2014 Kelly Kearney

% Set up

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

% Flag for nemuro-only (as opposed to water column ecosystem)

p.addParamValue('isnem',        true,                      @(x) islogical(x) && isscalar(x));

% Both

p.addParamValue('NemParam',     [],                        @isstruct);
p.addParamValue('bnem0',        []);
p.addParamValue('odesolver',    {'euler'});
p.addParamValue('reroute',      cell(0,5),                 @(x) isempty(x) || (iscell(x) && size(x,2)==5));

p.addParamValue('remineff',     0.5,                       @(x) isnumeric(x) && isscalar(x));
p.addParamValue('kliglo',       1.0e12,                    @(x) isnumeric(x) && isscalar(x));
p.addParamValue('klighi',       1.0e8,                     @(x) isnumeric(x) && isscalar(x));
p.addParamValue('iofescav',     10,                        @(x) isnumeric(x) && isscalar(x));
p.addParamValue('alphascav',    15.0/86400/365,            @(x) isnumeric(x) && isscalar(x));
p.addParamValue('felig2don',    0,                         @(x) isnumeric(x) && isscalar(x));
p.addParamValue('ligbkg',       1e-9,                      @(x) isnumeric(x) && isscalar(x));
p.addParamValue('kfe2n',        [3.0 6.0]*1e-6*106/16,     @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('Kfe',          [1e-10 5e-10],             @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('fe2nmax',      [50.0 500.0]*1e-6*106/16,  @(x) isnumeric(x) && isvector(x) && length(x)==2);
p.addParamValue('fe2nupfac',    15e-6,                     @(x) isnumeric(x) && isscalar(x));
p.addParamValue('scavthresh',   1e-6,                      @(x) isnumeric(x) && isscalar(x));
p.addParamValue('scavfac',      0,                         @(x) isnumeric(x) && isscalar(x));    

p.addParamValue('diapause',     false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('dday',         [120 134 243],             @(x) isnumeric(x) && isvector(x) && length(x)==3);
p.addParamValue('dfrac',        0.9,                       @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
% p.addParamValue('grazeatdepth', true,                      @(x) islogical(x) && isscalar(x));

% Wce-only

p.addParamValue('nomix',        false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('ecosimpp',     false,                     @(x) islogical(x) && isscalar(x));
p.addParamValue('Ewein',        [],                        @isstruct);
p.addParamValue('types',        []);
p.addParamValue('mld',          -50,                       @(x) isscalar(x) && x<=0);
p.addParamValue('temp',         0,                         @(x) isscalar(x));
p.addParamValue('kgra',         0.0693,                    @(x) isvector(x));
p.addParamValue('m0exp',        2,                         @(x) isnumeric(x) && (isscalar(x) || isvector(x)));
% p.addParamValue('predatdepth',  true,                      @(x) islogical(x) && isscalar(x));
p.addParamValue('preyvis',      [0 1; 5000 1],             @(x) isnumeric(x)); 

% Nemuro-only

p.addParamValue('ivlev',        'orig',                    @(x) ischar(x) && ismember(x, {'single', 'multi', 'orig', 'mishmash'}));
p.addParamValue('grmax',        [],                        @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 11));
p.addParamValue('thresh',       [],                        @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 11));
p.addParamValue('p',            [],                        @(x) isempty(x) || (isnumeric(x) && isequal(size(x), [11 11])));

% Parse

p.parse(varargin{:});
BioIn = p.Results;

% Check that required inputs are here

if BioIn.isnem
    nodefault = {'NemParam','bnem0'};
else
    nodefault = {'Ewein', 'NemParam', 'types', 'bnem0'};
end

tf = ismember(p.UsingDefaults, nodefault);
missing = sprintf('%s, ', p.UsingDefaults{tf});
if any(tf)
    error('Missing input for wce/nemuro module: %s', missing(1:end-2));
end

