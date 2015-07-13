function varargout = nemuroiron(action, varargin)
%NEMUROIRON NEMURO biological module with iron
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
%               (mol/l), where columns are as follows: 
%               1:  TPS, Small Phytoplankton
%               2:  TPL, Large Phytoplankton
%               3:  TZS, Small Zooplankton
%               4:  TZL, Large Zooplankton
%               5:  TZP, Pradatory Zooplankton
%               6:  TNO3, Nitrate
%               7:  TNH4, Ammmonium
%               8:  TPON, Particulate Organic Nitrogen
%               9:  TDON, dissolved Organic Nitrogen
%               10: TSiOH4, Silicate
%               11: TOpal, Particulate Opal
% 
%   bioinitz:   nz x 1 aaray, depths corresponding to rows of bioinit
% 
%   odesolver:
% 
%   reroute:    
% 
%   *nemparam*: scalar(s), any/all of the NEMURO parameters. See
%               nemuroinputparser.m for a full list of these parameters and
%               their default values.  Any parameter not included
%               explicitly as an input will be assigned the default value
%               according to the 'NEMURO Version 1.f90' set of parameters.

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

p = inputParser;

p.addParamValue('bioinit', nan(2,11), @(x) isequal(size(x,2), 11)); % Initial concentrations for biomass
p.addParamValue('bioinitz', [0 max(-Grd.z)]);
p.addParamValue('odesolver', {'ode4', 'ode45'});
p.addParamValue('reroute', cell(0,5), @(x) isempty(x) || (iscell(x) && size(x,2)==5));

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

In = mergestruct(p.Results, p.Unmatched);

[Biovars, Other] = nemuroinputparser(In);
if ischar(Other.odesolver)
    Biovars.odesolver = {Other.odesolver};
else
    Biovars.odesolver = Other.odesolver;
end

if isempty(Other.reroute)
    Biovars.reroute = cell(0,5);
else
    Biovars.reroute = Other.reroute;
end
    
defaultbio = [0.1 0.1 0.1 0.1 0.1 5.0 0.1 0.1 0.1 10 0.01] * 1e-6;

S = warning('off', 'MATLAB:interp1:NaNinY');
bio = interp1(In.bioinitz, In.bioinit, -Grd.z);
warning(S);

for ib = 1:size(bio,2)
    if all(isnan(bio(:,ib)))
        bio(:,ib) = defaultbio(ib);
    end
end

% All variables mixed

ismixed = true(1, size(bio,2));

% Names

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

%     {...
%     'Fe'        'Dissolved Iron'                'umolFe/l'
%     'PLsi'      'Large phytoplankton Silica'    'molSi/l'
%     'PSfe'      'Small phytoplankton Iron'      'umolFe/l'
%     'PLfe'      'Large phytoplankton Iron'      'umolFe/l'};

% No bottom forcing

bottomval = nan(size(ismixed));
% bottomval = nan(size(ismixed));
% bottomval([6 10]) = 0;

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

% Intermediate fluxes

[short, long] = findnemflux(nemuroflexinput(Biovars));

% Add rerouted fluxes

rerouteshort = cellfun(@(a,b,c) sprintf('%s_%02d_%02d', a,b,c), ...
    Biovars.reroute(:,1), Biovars.reroute(:,2), Biovars.reroute(:,4), 'uni', 0);
reroutelong = cellfun(@(a,b,c) sprintf('%s: %d to %d', a,b,c), ...
    Biovars.reroute(:,1), Biovars.reroute(:,2), Biovars.reroute(:,4), 'uni', 0);

short = [short; rerouteshort];
long  = [long;  reroutelong];

unit = repmat({'mol/l/s'}, size(short));

diagnames = [diagnames; [short long unit]];


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         

odefun = @(t,y,P) nemuroironode(t,y,Biovars,-z,dz,meanqi,temperature);

[newbio, db, Flx, Diag, badthings] = integratebio(odefun, t, dt, oldbio, [], Biovars.odesolver{:});

if any(badthings(:))
    error('Failed to integrate biology');
end

% Diagnostics

ndiag = 12;
if isempty(Diag) % If use any solver other than euler
    diag = zeros(length(z), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(2, diag{:});
end


nemflux = findnemflux(nemuroflexinput(Biovars), 'list');
nemflux{1} = [nemflux{1} Biovars.reroute(:,1)'];            % add rerouted
nemflux{2} = [nemflux{2} cell2mat(Biovars.reroute(:,2)')];  % add rerouted
nemflux{3} = [nemflux{3} cell2mat(Biovars.reroute(:,4)')];  % add rerouted

nfx = length(nemflux{2});
nz = size(newbio,1);
fluxes = zeros(nz, nfx);
for ifx = 1:nfx
    fluxes(:,ifx) = Flx(1).(nemflux{1}{ifx})(nemflux{2}(ifx), nemflux{3}(ifx),:);
end

diag = [diag fluxes];

%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% PON and Opal settle

biosettle = zeros(1,11);
biosettle(1,8) = -Biovars.setVPON; 
biosettle(1,11) = -Biovars.setVOpal;

wsink = ones(size(z)) * biosettle;
