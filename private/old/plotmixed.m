function h = plotmixed(file, varargin)
%PLOTMIXED Plots results of the mixed_layer model
%
% h = plotmixed(file, field1, field2, ...)
% h = plotmixed(file, splitdepth, ...)
%
% Input variables:
%
%   file:       mixed_layer output file
%
%   field:      output variables to be plotted, corresponding to a variable
%               names in the output file, or an arbitrary matrix of size
%               ndepth x ntime, ndepth+1 x ntime, or 1 x ntime.
%               
%
%   splitdepth: depth to separate surface and deep layers.  The surface
%               layers will be displayed in the upper quarter of the
%               figure; the bottom layers in the lower three quarters.
%
% Output variables:
%
%   h:          structure with handles to figures and axes

% Copyright 2008 Kelly Kearney

%------------------------------
% Check input
%------------------------------

if ~strcmp(file(end-2:end), '.nc')
    file = [file '.nc'];
end
if ~exist(file, 'file')
    error('Results file not found');
end

Info = nc_info(file);

% Some versions of nc_info use DataSet, some Dataset

fld = fieldnames(Info);
isds = strcmpi('dataset', fld);
dsfld = fld{isds};

% Extract depth and time variables

z = -nc_varget(file, 'depth');
zp = -nc_varget(file, 'depth_edge');
mtime = nc_varget(file, 'middate');
stime = nc_varget(file, 'startdate');
etime = nc_varget(file, 'enddate');

% Parse depth to split surface from deep

isname = cellfun(@(x) ischar(x) && ismember(x, {Info.(dsfld).Name}), varargin);
issplit = cellfun(@(x) isnumeric(x) && isscalar(x), varargin);
ismatrix = cellfun(@isnumeric, varargin) & ~issplit;

if any(~isname & ~issplit & ~ismatrix)
    error('Inputs must be variable names, split depth, or numeric arrays');
end

matidx = find(ismatrix);
newfield = arrayfun(@(x) sprintf('userdefined%d', x), 1:length(matidx), 'uni', 0);
    
if any(issplit)
    splitdepth = varargin{issplit};
else
    splitdepth = .25 * max(z);
end

vars = varargin;
vals = varargin;

vars(ismatrix) = newfield;

for iv = find(isname)
    vals{iv} = nc_varget(file, vars{iv})';
end

vars(issplit) = [];
vals(issplit) = [];

% Remove nulls

for iv = 1:length(vals)
    isnull = vals{iv} > 1e36;
    vals{iv}(isnull) = NaN;
end

%------------------------------
% Setup
%------------------------------

ndim = cellfun(@ndims, vals);

nvar = length(vars);

if nvar < 1
    error('No mixed_layer variable specified');
end

nz  = length(z);
nzp = length(zp);
nt  = length(mtime);

time = mtime;
time2 = [stime; etime(end)];

isshallow = z <= splitdepth;

mr = .15;
props = {'mr', mr, 'sv', .01};



%------------------------------
% Plot figures
%------------------------------

figcnt = 0;

for ivar = 1:nvar
    
    % Plot vector data as line plot
    
    if isequal(size(vals{ivar}), [1 nt])
        
        figcnt = figcnt + 1;
        
        h.fig(figcnt) = figure;
        
        h.ax(figcnt,1) = axes;
        plot(time, vals{ivar});
        xlabel('Time');
        ylabel(vars{ivar});
        title(vars{ivar});
        tlabel(h.ax(figcnt,1));
        
    % Plot 2-d arrays as pcolor plot
        
    elseif isequal(size(vals{ivar}), [nz nt])
        
        figcnt = figcnt + 1;
        
        h.fig(figcnt) = figure;
        
        h.ax(figcnt,1) = subaxis(3,1,1, props{:});
        pcolorpad(time2, zp, vals{ivar});
        shading flat;
        
        h.ax(figcnt,2) = subaxis(3,1,1,2,1,2, props{:});
        pcolorpad(time2, zp, vals{ivar});
        shading flat;
        
        set(h.ax(figcnt,1), 'ylim', [0 splitdepth]);
        set(h.ax(figcnt,2), 'ylim', [splitdepth max(zp)]);
        
        clims = minmax(vals{ivar}(~isinf(vals{ivar})));
%         clims = [0 clims(2)];
        labelpcolor(h.ax(figcnt,:), vars{ivar}, clims, mr);
        tlabel(h.ax(figcnt,:), 'whichaxes', 'last');
       
        
    elseif isequal(size(vals{ivar}), [nzp nt])
        
        figcnt = figcnt + 1;
        
        h.fig(figcnt) = figure;
        
        h.ax(figcnt,1) = subaxis(3,1,1, props{:});
        pcolorpad(time, zp(isshallow), vals{ivar}(isshallow,:));
        shading flat;
        
        h.ax(figcnt,2) = subaxis(3,1,1,2,1,2, props{:});
        pcolorpad(time, zp(~isshallow), vals{ivar}(~isshallow,:));
        shading flat;
        
        clims = minmax(vals{ivar}(~isinf(vals{ivar})));
%         clims = [0 clims(2)];
        labelpcolor(h.ax(figcnt,:), vars{ivar}, clims, mr);
        tlabel(h.ax(figcnt,:), 'whichaxes', 'last');
    
    end
end

% Switch to zbuffer renderer to correct some opengl weirdness

set(h.fig, 'renderer', 'zbuffer');

%------------------------------
% Pcolor plot with last row and
% column displayed
%------------------------------

function h = pcolorpad(x,y,z)

if isequal(size(z), [length(y) length(x)])
    x2 = [x(:); NaN];
    y2 = [y(:); NaN];
elseif isequal(size(z), [length(y)-1 length(x)-1])
    x2 = x;
    y2 = y;
end

[nrow, ncol] = size(z);
z2 = nan(nrow+1, ncol+1);
z2(1:nrow, 1:ncol) = z;

h = pcolor(x2, y2, z2);

%------------------------------
% Add axis labels and colorbar
% to both pcolor axes
%------------------------------

function labelpcolor(ax, ttl, clims, mr)

if diff(clims) == 0
    clims = [-1 1] + clims(1); % In case of no variation
elseif any(isnan(clims))
    clims = [-1 1]; % In case of all NaN
end

set(ax(1), 'xticklabel', '');
set(ax, 'ydir', 'reverse', ...
        'clim', clims, ...
        'layer', 'top');

pos = cell2mat(get(ax, 'position'));
yposlim = minmax([pos(:,2) pos(:,2)+pos(:,4)]);
xposmax = max(pos(:,1)+pos(:,3));

try
    cols = cptcmap('seminf-haxbycont'); % I like this colormap better than defaults, but most people won't have it
catch
    cols = jet(64);
end

colormap(cols);

cb = colorbar('peer', ax(1), 'position', [xposmax+mr*.25 yposlim(1) mr*.3 diff(yposlim)]);

[hy,ht,hax] = suplabel('axes', ax, 'ylabel', 'Depth (m)', ...
                                      'title',  ttl);                                 

set(hax, 'hittest', 'off');
