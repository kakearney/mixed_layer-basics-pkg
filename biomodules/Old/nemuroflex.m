function varargout = nemuroflex(action, varargin)
%NEMUROFLEX A slightly-more-flexible NEMURO biological module
%
% Experimental, was trying to do a rewrite of NEMURO while allowing more
% than three zooplankton groups.  

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

Biovars = nemuroflexinput(In);
nb = length(Biovars.alpha);

p = inputParser;

p.addParamValue('bioinit', nan(2,nb), @(x) isequal(size(x,2), nb)); % Initial concentrations for biomass
p.addParamValue('bioinitz', [0 max(-Grd.z)]);

p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(In);

In = p.Results;

defaultbio = [ones(1,nb)*1e-7 [5.0 0.1 0.1 0.1 10 0.01] * 1e-6];

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

nex = nb - 11;
if nex > 0
    exshort = cellstr(num2str((1:nex)', 'EX%d')); 
    exlong  = cellstr(num2str((1:nex)', 'Extra critter %d')); 
    [exunit{1:nex}] = deal('molN/l');
    names = [names(1:5,:); [exshort exlong exunit']; names(6:end,:)];
end
    
% No bottom forcing

bottomval = nan(size(ismixed));

% Diagnostics

% diagnames = {
%     'PSlightlim',   'Light limitation (small)',         'no units'
%     'PLlightlim',   'Light limitation (large)',         'no units'
%     'PSno3lim',     'Nitrate limitation (small)',       'no units'
%     'PLno3lim',     'Nitrate limitation (large)',       'no units'
%     'PSnh4lim',     'Ammonium limitation (small)',      'no units'
%     'PLnh4lim',     'Ammonium limitation (large)',      'no units'
%     'PStemplim',    'Temperature limitation (small)'    'no units'
%     'PLtemplim',    'Temperature limitation (large)'    'no units'
%     'PLsilim',      'Silica limitation (large)'         'no units'
%     'I',            'Irradiance'                        'W m^-2'
%     'kappa',        'Attenuation coefficient'           'm^-1'
%     'kp'            'Attenuation self-shading only',    'm^-1'};


% [nt,nb] = size(bio);
% bidx = (1:nb)';
% dshort = cellstr([num2str(bidx, 'qb%02d'); num2str(bidx, 'pb%02d')]); 
% dlong = cellstr([num2str(bidx, 'Q/B for %02d'); num2str(bidx, 'P/B for %02d')]);
% [dunit{1:nb*2}] = deal('s^-1');
% dunit = dunit';
% diagnames = [dshort dlong dunit];
diagnames = cell(0,3);

% Intermediate fluxes:

[fluxshort, fluxlong] = findnemflux(Biovars, 'flex');
[fluxunit{1:length(fluxshort)}] = deal('mol/l/s');

diagnames = [diagnames; [fluxshort fluxlong fluxunit']];

% nemflux = {...
%     'Gpp', [6 7 6 7 10], [1 1 2 2 12]
%     'Gra', [1 1 2 2 3 3 4], [3 4 4 5 4 5 5]
%     'Res', [1 2 1 2 12], [6 6 7 7 10]
%     'Exc', [1 2 3 4 5 12], [9 9 7 7 7 10]
%     'Ege', [3 4 5 12], [8 8 8 11]
%     'Mor', [1 2 3 4 5 12], [8 8 8 8 8 11]
%     'Dec', [8 9 8 11 7], [9 7 7 10 6]};
% 
% for ifx = 1:size(nemflux,1)
%     for ic = 1:length(nemflux{ifx,2})
%         short = sprintf('%s_%02d_%02d', nemflux{ifx,1}, nemflux{ifx,2}(ic), nemflux{ifx,3}(ic));
%         long = sprintf('%s: %d to %d', nemflux{ifx,1}, nemflux{ifx,2}(ic), nemflux{ifx,3}(ic));
%         diagnames = [diagnames; {short long 'mol/l/s'}];
%     end
% end


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                         
% Integrate using Runge-Kutta solver over the full timestep

[tout,newbio,db,Flx,Diag] = odewrap(...
    @ode4, @nemuroflexode, [t t+dt], oldbio, [], Biovars, -z, dz, ...
    meanqi, temperature);
                         
newbio = endonly(newbio);

% If something weird happens stepping over full timestep, use variable-step
% solver to try instead

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    [tout,newbio,db,Flx,Diag] = odewrap(...
        @ode45, @nemuroflexode, [t t+dt], oldbio, [], Biovars, -z, dz, ...
        meanqi, temperature);
    newbio = endonly(newbio);
end

% If still no good, error

if any(isnan(newbio(:)) | isinf(newbio(:)) | newbio(:) < 0)
    error('Biology problems');
end

% Diagnostics

diag = struct2cell(Diag(1)); % use beginning of time step
diag = cat(2, diag{:});
% diag = zeros(length(z),0);

nemflux = findnemflux(Biovars, 'list', 'flex');

nfx = length(nemflux{2});
nz = size(newbio,1);
fluxes = zeros(nz, nfx);
for ifx = 1:nfx
    fluxes(:,ifx) = Flx(1).(nemflux{1}{ifx})(nemflux{2}(ifx), nemflux{3}(ifx),:);
end

% nemflux = {...
%     'Gpp', [6 7 6 7 10], [1 1 2 2 12]
%     'Gra', [1 1 2 2 3 3 4], [3 4 4 5 4 5 5]
%     'Res', [1 2 1 2 12], [6 6 7 7 10]
%     'Exc', [1 2 3 4 5 12], [9 9 7 7 7 10]
%     'Ege', [3 4 5 12], [8 8 8 11]
%     'Mor', [1 2 3 4 5 12], [8 8 8 8 8 11]
%     'Dec', [8 9 8 11 7], [9 7 7 10 6]};
% 
% fluxes = zeros(size(diag,1), 37);
% count = 0;
% for ifx = 1:size(nemflux,1)
%     for ic = 1:length(nemflux{ifx,2})
%         count = count + 1;
%         fluxes(:,count) = Flx(1).(nemflux{ifx,1})(nemflux{ifx,2}(ic), nemflux{ifx,3}(ic),:);
%     end
% end

diag = [diag fluxes];

%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

% PON and Opal settle

% biosettle = zeros(1,11);
% biosettle(1,8) = -Biovars.setVPON; 
% biosettle(1,11) = -Biovars.setVOpal;

wsink = ones(size(z)) * -Biovars.settle';
