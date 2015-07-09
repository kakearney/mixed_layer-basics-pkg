function h = plotforcing(varargin)
%PLOTFORCING Plots heat, wind, and T/S forcing
%
% h = plotforcing
% h = plotforcing(In)
% h = plotforcing(wind_input, heat_input, ts_input)
%
% This function plots the forcing data used for a mixed_layer simulation.
% If no input is supplied, it plots the forcing used for a default
% simulation.
%
% Input variables:
%
%   wind_input: Wind forcing data.  This is an n x 8 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               east-west wind speed (i.e. u), north-south wind speed (i.e.
%               v). Speeds in m/s.  By default, a dateset representing 1976
%               observations on the Scotia Shelf are used.
%
%   heat_input: Heat forcing data.  This is an n x 9 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               incident solar radiation (Qi), air temperature, and dew
%               point temperature.  Radiation is in W/m^2 and all
%               temperatures are in deg C. By default, a dateset
%               representing 1976 observations on the Scotia Shelf are
%               used.
%
%  ts_input:    Initial temperature and salinity profiles for simulation.
%               Data is an n x 3 matrix with columns representing depth (m,
%               negative down), temperature (deg C), and salinity (psu).
%               By default, a dateset representing 1976 observations on the
%               Scotia Shelf are used.
%
%   In:         structure with fields 'wind_input', 'heat_input', and
%               'ts_input', holding arrays as described above
%
% Output variables
%
%   h:          structure of handles

% Copyright 2009 Kelly Kearney

%---------------------------
% Check input
%---------------------------

if nargin == 3
    In.wind_input = varargin{1};
    In.heat_input = varargin{2};
    In.ts_input   = varargin{3};
elseif nargin == 0
    In = parseinput('blah');
elseif nargin == 1
    In = varargin{1};
end 

%---------------------------
% Plot
%---------------------------

h.fig = figure('color', 'w');
h.panel(1) = axes('position', [.1 .1 .55 .35]);  % wind
h.panel(2) = axes('position', [.75 .3 .2 .6]);  % ts
h.panel(3) = axes('position', [.1 .5 .55 .4]);  % heat

% Plot wind 

axes(h.panel(1));
[h.lines(1:2), h.axes(1:2)] = plotses(datenum(In.wind_input(:,1:6)), In.wind_input(:,7:8));

try
    tlabel(h.axes(1), 'keeplimits');%, 'FixLow', 4);
catch
    tlabelold(h.axes(1), 'keeplimits');
end

% Plot heat forcing

axes(h.panel(3));
[h.lines(3:5), h.axes(3:5)] = plotses(datenum(In.heat_input(:,1:6)), In.heat_input(:,7:9));

% try
    tlabel(h.axes(3), 'keeplimits');%, 'FixLow', 4);
% catch
%     tlabelold(h.axes(3), 'keeplimits')
% end

% Plot temperature and salinity

z = In.ts_input(:,1);
y = In.ts_input(:,2:3);

axes(h.panel(2));
[h.lines(6:7), h.axes(6:7)] = plots(z,y, 'top');

% Colors and legend

try
    cols = cptcmap('Dark2_08', 'ncol', 8);
catch
    cols = hsv(7);
end
    
labels = {...
    'U wind (m/s)' 
    'V wind (m/s)' 
    'Radiation (W/m^2)' 
    'Air temp (deg C)' 
    'Dew point temp. (deg C)' 
    'Temperature (deg C)'
    'Salinity (ppt)'};

for il = 1:7
    set(h.lines(il), 'color', cols(il,:));
end

h.legend = legend(h.lines, labels);
set(h.legend, 'fontsize', 10, 'position', [.75 .1 .2 .2]);
legend(h.legend, 'boxoff');