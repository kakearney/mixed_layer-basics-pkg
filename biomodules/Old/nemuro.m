function varargout = nemuro(action, varargin)
%NEMURO NEMURO biological module
%
% See biomodule.m for full syntax details
%
% This module runs a version of the NEMURO model, based on the NEMURO
% Version 1 code available from PICES.  It includes 11 biological state
% variables (2 phytoplankton, 3 zooplankton, and N- and Si-cycles).
%
% User-specified input variables (passed to mixed_layer as parameter/value
% pairs)
%
%   bioinit:    nz x 11 array, initial concentrations of state variables
%               (molN/l), where columns are as follows: 
%               1:  PS, Small Phytoplankton
%               2:  PL, Large Phytoplankton
%               3:  ZS, Small Zooplankton
%               4:  ZL, Large Zooplankton
%               5:  ZP, Pradatory Zooplankton
%               6:  NO3, Nitrate
%               7:  NH4, Ammmonium
%               8:  PON, Particulate Organic Nitrogen
%               9:  DON, dissolved Organic Nitrogen
%               10: SiOH4, Silicate
%               11: Opal, Particulate Opal
% 
%   bioinitz:   nz x 1 array, depths corresponding to rows of bioinit
% 
%   odesolver:  cell array of strings, indicating which solvers to use.  If
%               the first one fails to integrate a timestep (i.e. causes
%               something to become negative, NaN, or Inf), the next one is
%               tried.  See integratebio.m for choices. [{'ode4', 'ode45'}]
% 
%   reroute:    n x 5 cell array.  This allows you to reroute fluxes from
%               the original path defined in the nemuro model.  Each row
%               desribes as change, with column as follows:
%               col 1:  name of flux (Gpp, Gra, Res, Exc, Ege, Mor, Dec,
%                       Nit, or Ufe)
%               col 2:  index of original source group
%               col 3:  index of original sink group
%               col 4:  index of new sink group
%               col 5:  fraction of flux to reroute
% 
%   *nemparam*: scalar(s), any/all of the NEMURO parameters. See
%               nemuroinputparser.m for a full list of these parameters and
%               their default values.  Any parameter not included
%               explicitly as an input will be assigned the default value
%               according to the 'NEMURO Version 1.f90' set of parameters.
%
%   Fe:         1 x 1 structure of iron-related parameters.  If included,
%               iron limitation based on Fiechter et al., 2009 is added to
%               the model, and four additional biological state variables
%               are added:
%               12: Fe, Dissolved Iron (umolFe/l)
%               13: PLsi, Large phytoplankton Silica (molSi/l)
%               14: PSfe, Small phytoplankton Iron (umolFe/l)
%               15: PLfe, Large phytoplankton Iron (umolFe/l)
%
%               The structure holds the following fields:
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
%                       concentration (umol/l)

% Copyright 2008-2011 Kelly Kearney

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

%----------------------
% Parse input
%----------------------

p = inputParser;

p.addParamValue('bioinit', nan(2,11), @(x) isequal(size(x,2), 11)); % Initial concentrations for biomass
p.addParamValue('bioinitz', [0 max(-Grd.z)]);
p.addParamValue('odesolver', {'ode4', 'ode45'});
p.addParamValue('reroute', cell(0,5), @(x) isempty(x) || (iscell(x) && size(x,2)==5));

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

In = mergestruct(p.Results, p.Unmatched);

%----------------------
% Reformat/check input
%----------------------

% Separate NEMURO paramters from other stuff

[Biovars, Other] = nemuroinputparser(In);

% Set ODE solver

if ischar(Other.odesolver)
    Biovars.odesolver = {Other.odesolver};
else
    Biovars.odesolver = Other.odesolver;
end

% Look for rerouted fluxes

if isempty(Other.reroute)
    Biovars.reroute = cell(0,5);
else
    Biovars.reroute = Other.reroute;
end

% Check whether iron should be included

if isfield(Other, 'Fe')
    useiron = true;
    Biovars.Fe = Other.Fe;
    
%     flds = {'tfe', 'a', 'b', 'kfe', 'frem', 'cpern'};
    flds = {'kfe', 'frem'};
    for ifld = 1:length(flds)
        if isscalar(Biovars.Fe.(flds{ifld}))
            Biovars.Fe.(flds{ifld}) = [1 1] * Biovars.Fe.(flds{ifld});
        end
    end
else
    useiron = false;
end

Biovars.useiron = useiron;

    
%----------------------
% Initial biomass
%----------------------

defaultbio = [0.1 0.1 0.1 0.1 0.1 5.0 0.1 0.1 0.1 10 0.01] * 1e-6;

S = warning('off', 'MATLAB:interp1:NaNinY');
bio = interp1(In.bioinitz, In.bioinit, -Grd.z);
warning(S);

for ib = 1:size(bio,2)
    if all(isnan(bio(:,ib)))
        bio(:,ib) = defaultbio(ib);
    end
end

if useiron
    bio(:,12) = interp1(Biovars.Fe.fe0(:,1), Biovars.Fe.fe0(:,2), -Grd.z);
    bio(:,13) = bio(:,2) .* Biovars.RSiN;
    
    cperfe = 33e3;
    cpern = 6.625;
    fepern = cpern./cperfe; %  mol Fe / mol N
    fepern = fepern * 1e6;  % umol Fe / mol N
    
    bio(:,14:15) = bio(:,1:2) * fepern; % umol Fe/l
end

%----------------------
% Mixing and forcing
%----------------------

% All variables mixed

ismixed = true(1, size(bio,2));

% No bottom forcing

bottomval = nan(size(ismixed));

%----------------------
% Variable names
%----------------------

% State variables

names = {
    'PS',       'Small Phytoplankton',          'molN/l'
    'PL',       'Large Phytoplankton',          'molN/l'
    'ZS',       'Small Zooplankton',            'molN/l'
    'ZL',       'Large Zooplankton',            'molN/l'
    'ZP',       'Pradatory Zooplankton',        'molN/l'
    'NO3',      'Nitrate',                      'molN/l'
    'NH4',      'Ammmonium',                    'molN/l'
    'PON',      'Particulate Organic Nitrogen', 'molN/l'
    'DON',      'dissolved Organic Nitrogen',   'molN/l'
    'SiOH4',    'Silicate',                     'molSi/l'
    'Opal',     'Particulate Opal',             'molSi/l'};

if useiron
    names = [names; {...
    'Fe'        'Dissolved Iron'                'umolFe/l'
    'PLsi'      'Large phytoplankton Silica'    'molSi/l'
    'PSfe'      'Small phytoplankton Iron'      'umolFe/l'
    'PLfe'      'Large phytoplankton Iron'      'umolFe/l'}];
end
    
% Diagnostics

diagnames = {
    'PSlightlim',   'Light limitation (small)',         'no units'
    'PLlightlim',   'Light limitation (large)',         'no units'
    'PSno3lim',     'Nitrate limitation (small)',       'no units'
    'PLno3lim',     'Nitrate limitation (large)',       'no units'
    'PSnh4lim',     'Ammonium limitation (small)',      'no units'
    'PLnh4lim',     'Ammonium limitation (large)',      'no units'
    'PStemplim',    'Temperature limitation (small)'    'no units'
    'PLtemplim',    'Temperature limitation (large)'    'no units'
    'PLsilim',      'Silica limitation (large)'         'no units'
    'I',            'Irradiance'                        'W m^-2'
    'kappa',        'Attenuation coefficient'           'm^-1'
    'kp'            'Attenuation self-shading only',    'm^-1'};

if useiron
    diagnames = [diagnames; {...
    'PSfelim'       'Iron limitation (small)'           'no units'
    'PLfelim'       'Iron limitation (large)'           'no units'}];
%     'PSR0'          'Empirical Fe:C ratio (small)'      'umolFe/molC'
%     'PSR'           'Realized Fe:C ratio (small)'       'umolFe/molC'
%     'PLR0'          'Empirical Fe:C ratio (large)'      'umolFe/molC'
%     'PLR'           'Realized Fe:C ratio (large)'       'umolFe/molC'}];
end
   

% Intermediate fluxes

if useiron
    fluxlist = fluxindices('nemuroiron');
else
    fluxlist = fluxindices('nemuro');
end

fluxlist = [fluxlist; Biovars.reroute(:,[1 2 4])];
nd = size(fluxlist,1);

fluxnames = cell(nd,3);
for id = 1:nd
    fluxnames{id,1} = sprintf('%s_%02d_%02d', fluxlist{id,:});
    fluxnames{id,2} = sprintf('%s: %02d to %02d', fluxlist{id,:});
    fluxnames{id,3} = 'mol/l/s (or umol/l/s for Fe)'; % TODO fix this so gives proper N,Si,or Fe unit
end

Biovars.nemflux = fluxlist;

% [short, long] = findnemflux(nemuroflexinput(Biovars));

% Add rerouted fluxes

% rerouteshort = cellfun(@(a,b,c) sprintf('%s_%02d_%02d', a,b,c), ...
%     Biovars.reroute(:,1), Biovars.reroute(:,2), Biovars.reroute(:,4), 'uni', 0);
% reroutelong = cellfun(@(a,b,c) sprintf('%s: %d to %d', a,b,c), ...
%     Biovars.reroute(:,1), Biovars.reroute(:,2), Biovars.reroute(:,4), 'uni', 0);
% 
% short = [short; rerouteshort];
% long  = [long;  reroutelong];

% unit = repmat({'mol/l/s'}, size(short));

diagnames = [diagnames; fluxnames];


% Lag times for copepod grazing

if isfield(Other, 'lag')
    nlag = ceil(Other.lag./In.dt);
end
Biovars.laggedbio = bio;


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
if Biovars.useiron
    odefun = @(t,y,P) nemuroironode(t,y,Biovars,-z,dz,meanqi,temperature);
else
    odefun = @(t,y,P) nemuroode(t,y,Biovars,-z,dz,meanqi,temperature);
end

[newbio, db, Flx, Diag, badthings] = integratebio(odefun, t, dt, oldbio, [], Biovars.odesolver{:});


if any(badthings(:))
    
    % Iron is the only thing that might go negative for non-numerical
    % reasons
    
%     if Biovars.useiron && ~any(any(badthings(:,[1:11 13:15])))
%         isneg = badthings(:,12);
%         fe = oldbio(:,12);
%         db(isneg) = fe(isneg)./dt;
%         newbio(isneg) = 
%         blah
%         
%     else
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

% Diagnostics

if Biovars.useiron
    ndiag = 18;
else
    ndiag = 12;
end
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

ng = size(oldbio,2);

% PON and Opal settle

biosettle = zeros(1,ng);
biosettle(1,8) = -Biovars.setVPON; 
biosettle(1,11) = -Biovars.setVOpal;

wsink = ones(size(z)) * biosettle;
