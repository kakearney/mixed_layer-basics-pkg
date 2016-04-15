function varargout = nemuro_nosourcesink(action, varargin)
%NEMURO_NOSOURCESINK Template biological module for mixed layer model
%
% For debugging, probably not stable anymore

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

Biovars = nemuroinputparser(In);

p = inputParser;

p.addParamValue('bioinit', nan(2,11), @(x) isequal(size(x,2), 11)); % Initial concentrations for biomass
p.addParamValue('bioinitz', [0 max(-Grd.z)]);

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

In = p.Results;

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
    'SiOH4',    'Silicate',                     'molN/l'
    'Opal',     'Particulate Opal',             'molN/l'};

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

% diagnames = cell(0,3);


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);

%***Debugging***                         
newbio = oldbio;
diag = zeros(length(z),12);
return
%***End Debugging***   

% Integrate using Runge-Kutta solver over the full timestep

[tout,newbio,db,Flx,Diag] = odewrap(@ode4, @nemuroode, [t t+dt], oldbio, [], ...
                        Biovars, -z, dz, meanqi, temperature);
                         
[newbio, diag{1}] = endonly(newbio, db);


% diag = [diag; struct2cell(Flx(end))]';

if t <= 43200
    save(sprintf('testfluxn_%05d', t));
end

% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio] = odewrap(@ode45, @nemuroode, [t t+dt], oldbio, [], ...
                            Biovars, -z, dz, meanqi, temperature);
    newbio = endonly(newbio);
end

% If still no good, error

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end

% Diagnostics

diag = struct2cell(Diag(1));    % Look at limiters from beginning of time step
diag = cat(2, diag{:});


%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% PON and Opal settle

biosettle = zeros(1,11);
biosettle(1,8) = -Biovars.setVPON; 
biosettle(1,11) = -Biovars.setVOpal;

wsink = ones(size(z)) * biosettle;
