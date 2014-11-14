function In = parseinput(varargin)   
%PARSEINPUT Parse and check all input variables
%
% In = parseinput(varargin) 
%
% This function checks all input variables passed to the mixed_layer
% model. It supplies default values where applicable, and runs some simple
% validation checks on the values of each variable. See mixed_layer for
% input variable list.
% 
% Input variables:
%
%	see mixed_layer
%
% Output variables:
%
% 	In:	structure holding all input variables, with default values substituted
% 		where applicable   

% Copyright 2012 Kelly Kearney 		
                                        
p = inputParser;

isnumscal = @(x) isnumeric(x) && isscalar(x);

% No partial matches; it can cause a seg fault w/ some of my options

if isprop(p, 'PartialMatching') % only R2013+
    ip.PartialMatching = false;
end

% The only required input: the output file name

p.addRequired('outputfile', @ischar);

% Model parameters and their defaults

p.addParamValue('dz', 5, @(x) isnumeric(x) && isvector(x) && size(x,2)==1);
p.addParamValue('zbot', -150, @(x) isnumeric(x) && isscalar(x) && x < 0);
p.addParamValue('syear', 1976, @(x) isscalar(x) || isequal(size(x), [1 6]));
p.addParamValue('eyear', 1976, @(x) isscalar(x) || isequal(size(x), [1 6]));
p.addParamValue('dt', 3*3600, isnumscal);
p.addParamValue('tarch', 24*3600, @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParamValue('whgt', 10, isnumscal);
p.addParamValue('Lat', 45, isnumscal);
p.addParamValue('krad1', 0.15, isnumscal);
p.addParamValue('prad1', 0.45, isnumscal);
p.addParamValue('krad2', 1.67, isnumscal);
p.addParamValue('alb', 0.079, isnumscal);
p.addParamValue('pgx', 1e-5, isnumscal);
p.addParamValue('pgy', 0, isnumscal);
p.addParamValue('kmol', 1e-4, isnumscal);
p.addParamValue('velocity_dissipation', 1/(3*86400), isnumscal);
p.addParamValue('srelaxtime', 3600*24*30, isnumscal);
p.addParamValue('trelaxtime', 3600*24*30, isnumscal);
p.addParamValue('tempfilesz', NaN, isnumscal);
p.addParamValue('brelaxtime', 3600*24*30, @isvector);
p.addParamValue('openbottom', false, @(x) isscalar(x) && islogical(x));
p.addParamValue('beginarchive', NaN, @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParamValue('endarchive',   NaN, @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParamValue('outputextension', [], @(x) isempty(x) || (iscell(x) && all(cellfun(@ischar,x)))); 
p.addParamValue('tempdir', [], @(x) isempty(x) || ischar(x));
p.addParamValue('cleanup', true, @(x) isscalar(x) && islogical(x));
p.addParamValue('verbose', true, @(x) isscalar(x) && islogical(x));
p.addParamValue('hotstartdn', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
p.addParamValue('hotstartload', '',@(x) ischar(x));
p.addParamValue('hotstartsave', 'mlhotstart', @(x) ischar(x));
p.addParamValue('iens', 1, isnumscal);
% p.addParamValue('nens', 1, isnumscal);
% p.addParamValue('newfile', true, @(x) isscalar(x) && islogical(x));

% Optional model parameters (default to empty)

p.addParamValue('biofun',  [], @(x) isempty(x) || isa(x, 'function_handle'));
p.addParamValue('srelax',  [], @isnumeric);
p.addParamValue('trelax',  [], @isnumeric);
p.addParamValue('tracerw', [], @isnumeric);
p.addParamValue('advheatfun', [], @(x) isempty(x) || isa(x, 'function_handle'));

% Placeholders for forcing datasets

p.addParamValue('wind_input', [], @(x) size(x,2) == 8);
p.addParamValue('heat_input', [], @(x) size(x,2) == 9);
p.addParamValue('ts_input', [], @(x) size(x,2) == 3 && all(x(:,1)<=0));

% Allow structure input

p.StructExpand = true;

% Allow unmatched inputs (for biological module inputs, which could be
% anything)

p.KeepUnmatched = true;

% Parse input

p.parse(varargin{:});
In = mergestruct(p.Results, p.Unmatched);

% Add default datasets if necessary

mlname = mfilename('fullpath');
mlpath = fileparts(fileparts(mlname));
defaultdir = fullfile(mlpath, 'defaultdata');

datasets = {'wind_input', 'heat_input', 'ts_input'};
isdefault = ismember(datasets, p.UsingDefaults);
if ~all(isdefault) && any(isdefault)
    missingdata = sprintf('%s, ', datasets{isdefault});
    missingdata = missingdata(1:end-2);
    warning('ML:missingdata', 'Missing datasets: %s.\nYour model forcing will be a mix a your data and default data\n(which will probably lead to some weird results)', missingdata);
end

if isdefault(1)
    B = load(fullfile(defaultdir, 'wind_input.mat'));  % loads n x 8 matrix of wind data, see above
    In.wind_input = B.wind_input;
end

if isdefault(2)
    B = load(fullfile(defaultdir, 'heat_input.mat')); % loads n x 9 matrix of heat input, see above
    In.heat_input = B.heat_input;
end

if isdefault(3)
    B = load(fullfile(defaultdir, 'ts_input.mat'));   % loads n x 3 matrix of temperature and salinity, see above
    In.ts_input = B.ts_input;
end

% Allow .nc extension to be left off of output file name

[blah,blah,ext] = fileparts(In.outputfile);
if isempty(ext) || ~strcmp(ext, 'nc')
    In.outputfile = [In.outputfile '.nc'];
end

% Indicators of optional parameters

In.hasbio     = ~isempty(In.biofun);
In.hassrelax  = ~isempty(In.srelax);
In.hastrelax  = ~isempty(In.trelax);
In.hasw       = ~isempty(In.tracerw);
In.hasadvheat = ~isempty(In.advheatfun);

% Check for agreement between archiving variables

In.tarch        = In.tarch(:);
In.beginarchive = In.beginarchive(:);
In.endarchive   = In.endarchive(:);

len = cellfun(@length, {In.tarch, In.beginarchive, In.endarchive});
maxlen = max(len);
if ~all(len == maxlen | len == 1)
    error('Archiving variables must agree in number');
end
if maxlen > 1
    if len(1) == 1
        In.tarch = ones(maxlen,1)*In.tarch;
    end
    if len(2) == 1
        In.beginarchive = ones(maxlen,1)*In.beginarchive;
    end
    if len(3) == 1
        In.endarchive = ones(maxlen,1)*In.endarchive;
    end
end

if isempty(In.outputextension)
    In.outputextension = cellstr(num2str((1:length(In.tarch))'));
else
    if length(In.outputextension) ~= length(In.tarch)
        error('Incorrect number of output extensions');
    end
end

% Hot start parameters (save or use)

if ~isempty(In.hotstartload) && ~exist(In.hotstartload, 'file')
    error('Hot start file not found');
end
    
    
    