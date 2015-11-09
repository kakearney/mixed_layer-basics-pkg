%% Examples for the mixed_layer model

%% Setting up your path
% 
% The mixed_layer code is spread over several folders (mostly because that
% makes it easier to maintain on my end... sorry about that).  To get the
% basics running, you need to add a few folders to your path:  
% 
%   addpath('mixed_layer-basics-pkg/mixed_layer');
%   addpath('mixed_layer-basics-pkg/mixed_layer/biomodules');
%   addpath('mixed_layer-basics-pkg/mixed_layer/biomodules_subs');
%   addpath('mixed_layer-basics-pkg/seawater_ver3_2');
%   addpath('mixed_layer-basics-pkg/mergestruct');
%
% (Alter the paths as necessary to wherever you place these folders).
%
% If you only plan to run the |cobalt_fweb| module, or if you're adding
% your own biology module, then that's all you need.  To play with the |np|
% or |npz| modules, add three more (although these are simple models, I
% developed them when I was still experimenting with ODE solvers, hence
% the extra dependencies):
%
%   addpath('mixed_layer-basics-pkg/odefixed');
%   addpath('mixed_layer-basics-pkg/odewrap');
%   addpath('mixed_layer-basics-pkg/photosynthesis');
%
% Most users (i.e. probably everyone but me) can ignore the rest of the
% folders, which are related to the wce (water 
% column ecosystem) and the nemurokak (which is the lower-trophics-only
% version of wce) modules.  The parameter setup for these modules can get
% very complicated, enough so that I have entire toolboxes worth of
% functions dedicated to the task for just my one region of interest.  I
% would really appreciate if you would email me if you intend to start
% using these options, so we can discuss potential collaborations (and I
% can give more specific suggestions on how to use this portion of the
% code).  That said, if you really want to give it a try, add the rest of
% the folders to your path:
%
%   addpath('mixed_layer-basics-pkg/ConsoleProgressBar');
%   addpath('mixed_layer-basics-pkg/cellstr2');
%   addpath('mixed_layer-basics-pkg/cprintf');
%   addpath('mixed_layer-basics-pkg/dirfull');
%   addpath('mixed_layer-basics-pkg/ecopathensemble');
%   addpath('mixed_layer-basics-pkg/ecopathlite');
%   addpath('mixed_layer-basics-pkg/legendflex');
%   addpath('mixed_layer-basics-pkg/parsepv');
%   addpath('mixed_layer-basics-pkg/regexpfound');
%   addpath('mixed_layer-basics-pkg/setgetpos_V1.2');



%% Running physics-only simulations
%
% The simplest way to run the |mixed_layer| model is to do so without any
% of the biological modules turned on.  The code comes with some default
% surface forcing datasets (1976 data for the Scotian Shelf), which can be
% used for some quick tests.  To run with all default variables, simply
% specify an output folder name.  

mixed_layer('default');

%%
% You might notice the warning messages under the Initializing step;
% they're letting us know that our forcing dataset doesn't quite span the
% whole simulation timeperiod, so it had to be extrapolated.  For tiny gaps
% like this, extrapolation is fine, but if you see gaps of more than a few
% hours, your input probably needs to be adjusted.  
%
% When run in single-simulation mode like this, the model creates two
% output files: one that holds the dimension variables (dimensions.nc), and
% one for everything else (sim0001.nc).  To read in data in the following
% examples, I use a little wrapper function called |ncreads.m|; you can
% find that <https://github.com/kakearney/ncreads-pkg here>.

if isempty(which('ncreads'))
    error('To run this script as is, you need to download ncreads: https://github.com/kakearney/ncreads-pkg');
end

Dim1  = ncreads(fullfile('default','dimensions.nc'));
Data1 = ncreads(fullfile('default','sim0001.nc'), 'temp', 'sal');

figure;
ax(1) = subplot(2,1,1);
pcolor(Dim1.middate, Dim1.depth, Data1.temp);
datetick('x', 'mmm');
shading flat;
colorbar;
title('Temperature');
ax(2) = subplot(2,1,2);
pcolor(Dim1.middate, Dim1.depth, Data1.sal);
datetick('x', 'mmm');
shading flat;
colorbar;
title('Salinity');
set(ax, 'layer', 'top');
set(ax(2), 'clim', [30 35]);

%%
% In practice, you'll want to use the appropriate surface forcings for
% your location and time period of interest.  Also, in almost all cases,
% you'll want to supply a salinity relaxation dataset; otherwise the water
% column salinity slowly homogenizes due to lack of freshwater inputs.
% This folder includes some example data for a more realistic run on the
% Eastern Scotian Shelf.

load heat_input_ess_1984_2007.mat;
load wind_input_ess_1984_2008.mat;
load ts_init_ess.mat;
load salt_clim_ess_1984_2008_smoothed;

mixed_layer('ess', ...
    'dt', 1800, ...
    'dz', 5, ...
    'zbot', -150, ...
    'syear', 1984, ...
    'eyear', 1985, ...
    'Lat', 45, ...
    'ts_input', ts_init_ess, ...
    'wind_input', wind_input_ess_1984_2008, ...
    'heat_input', heat_input_ess_1984_2007, ...
    'srelax', salt_clim_ess_1984_2008_smoothed, ...
    'srelaxtime', 259200);

%%

Dim2  = ncreads(fullfile('ess','dimensions.nc'));
Data2 = ncreads(fullfile('ess','sim0001.nc'), 'temp', 'sal');

figure;
ax(1) = subplot(2,1,1);
pcolor(Dim2.middate, Dim2.depth, Data2.temp);
datetick('x', 'mmmyy');
shading flat;
colorbar;
title('Temperature');
ax(2) = subplot(2,1,2);
pcolor(Dim2.middate, Dim2.depth, Data2.sal);
datetick('x', 'mmmyy');
shading flat;
colorbar;
title('Salinity');
set(ax, 'layer', 'top');

%% Adding biology
%
% In order to add biology, we need to turn on one of the biological modules
% via the |biofun| input, and supply any additional input parameters needed
% by the module.  See the function help for each module to figure out what
% these inputs are.
%
% Charlie's cobalt_fweb module is one of the easier ones to get started
% with, since most of its parameters are hard-coded.  It requires only one
% additional input, holding the initial concentration profiles for all the
% state variables; an example dataset is included in the examples folder.

load cobalt_init_ess.mat;

mixed_layer('ess_cobalt', ...
    'dt', 1800, ...
    'dz', 5, ...
    'zbot', -150, ...
    'syear', 1984, ...
    'eyear', 1985, ...
    'Lat', 45, ...
    'ts_input', ts_init_ess, ...
    'wind_input', wind_input_ess_1984_2008, ...
    'heat_input', heat_input_ess_1984_2007, ...
    'srelax', salt_clim_ess_1984_2008_smoothed, ...
    'srelaxtime', 259200, ...
    'biofun',@cobalt_fweb, ...
    'bio_input',cobalt_init_ess);

%%

Dim3  = ncreads(fullfile('ess_cobalt','dimensions.nc'));
vars = {'sp', 'lp', 'sz','mz','lz'};
nv = length(vars);
Data3 = ncreads(fullfile('ess_cobalt','sim0001.nc'), vars{:});

vpos = linspace(0.05, 0.95, nv+1);

figure;
for ii = 1:nv
    ax(ii) = axes('position', [0.1 vpos(ii), 0.8 (0.9/nv)*0.8]);
    pcolor(Dim3.middate, Dim3.depth, Data3.(vars{ii}));
    shading flat;
    title(vars{ii});
    colorbar;
end
set(ax(2:end), 'xticklabel', '');
datetick(ax(1), 'x', 'mmmyy', 'keeplimits');
set(ax, 'layer', 'top');

%% 
% In the above examples, we've been reusing a lot of the same inputs.  To
% make it easier to do this without constantly retyping, we can store the
% input variables in a sructure array, and pass that as input. 

In.dt = 1800;
In.dz = 5;
In.zbot = -150;
In.syear = 1984;
In.eyear = 1985;
In.Lat = 45;
In.ts_input = ts_init_ess;
In.wind_input = wind_input_ess_1984_2008;
In.heat_input = heat_input_ess_1984_2007;
In.srelax = salt_clim_ess_1984_2008_smoothed;
In.srelaxtime = 259200;

%%
% You can combine structure and parameter/value input.  This example reuses
% all the physical parameters listed above, but then runs the NP toy model
% instead.  Also, notice that if a parameter is repeated (as in the |eyear|
% parameter here, which is set both in the structure array and as a
% parameter/value pair), the last listed value is used.

mixed_layer('ess_np', In, 'n', cobalt_init_ess(:,[1 13]), ...
                          'p', cobalt_init_ess(:,[1 3]), ...
                          'kn', 0.1, ...
                          'loss', 0.05/86400, ...
                          'remin', 0.5, ...
                          'eyear', 1984, ...
                          'biofun', @np);

%%
Dim4  = ncreads(fullfile('ess_np','dimensions.nc'));
Data4 = ncreads(fullfile('ess_np','sim0001.nc'), 'N','P');

figure;
clear ax;
ax(1) = subplot(2,1,1);
pcolor(Dim4.middate, Dim4.depth, Data4.N);
datetick('x', 'mmmyy');
shading flat;
colorbar;
title('Nutrient');
ax(2) = subplot(2,1,2);
pcolor(Dim4.middate, Dim4.depth, Data4.P);
datetick('x', 'mmmyy');
shading flat;
colorbar;
title('Phytoplankton');
set(ax, 'layer', 'top');

%%
% In case you're wondering about the crazy combo of ODE solvers used in the
% NP and NPZ modules, that's mostly the result of some early
% experimentation when I was preparing to build a more complicated model...
% and that's also why it runs a bit slowly for such a simple model.

%% Ensemble simulations
%
% The mixed_layer model comes with the ability to run ensemble simulations,
% where the same basic model grid setup is used but the forcings or
% parameters are varied between simulations.  The |runmixedlayer.m|
% function can be used to make file management simpler in these cases.  For
% example, let's run a sensitivity study where we vary the half-saturation
% constant in the above example:

nens = 5;
kn = linspace(0.05, 3.5, nens);

Ens = In;
Ens.n = cobalt_init_ess(:,[1 13]);
Ens.p = cobalt_init_ess(:,[1 3]);
Ens.loss = 0.05/86400;
Ens.remin = 0.5;
Ens.eyear = 1984;
Ens.biofun = @np;

Ens = repmat(Ens, nens, 1);
for ii = 1:nens
    Ens(ii).kn = kn(ii);
end

runmixedlayer(Ens, 'name', 'ess_npsensitivity', 'verbose', false);

%%
% The output folder now contains simXXXX.nc files for the 5 different
% variations.

dir('ess_npsensitivity')

%%
% We could loop over all these files and read in the data.  But with larger
% ensembles and multiple variables, that can become a pain.  The
% postprocess option in runmixedlayer reformats the output data so that
% instead of one file per simulation with all the variables in it, you get
% one file per variable with all simulations.  The postprocessing options
% relies on the NCO utilities (http://nco.sourceforge.net/), so you'll need
% to have those installed and accessible via Matlab to use this option.

ncdisp(fullfile('ess_npsensitivity', 'sim0001.nc'));
%%
runmixedlayer(Ens, 'name', 'ess_npsensitivity', 'postprocess', true);

dir('ess_npsensitivity')
%%
ncdisp(fullfile('ess_npsensitivity', 'P.nc'));

%%
% Note that unlike mixed_layer itself, runmixedlayer never attempts to
% overwrite files.  If it finds a file already exists it skips to the next
% ensemble member, or to the postprocessing steps.  It also deals with
% errors differently; if one simulation crashes, the code catches the error
% and writes it to a log file, then moves on to the next step, rather than
% halting all simulations.
%
% We can now inspect the effects of changing $K_N$ on our NP model:

Dim5 = ncreads(fullfile('ess_npsensitivity','dimensions.nc'));
p = ncread(fullfile('ess_npsensitivity', 'P.nc'), 'P');

vpos = linspace(0.05, 0.95, nens+1);

figure;
for ii = 1:nens
    ax(ii) = axes('position', [0.1 vpos(ii), 0.8 (0.9/nens)*0.8]);
    pcolor(Dim5.middate, Dim5.depth, p(:,:,ii));
    shading flat;
    title(sprintf('K_N = %.3f', kn(ii)));
end
set(ax(2:end), 'xticklabel', '');
datetick(ax(1), 'x', 'mmmyy', 'keeplimits');
set(ax, 'layer', 'top', 'clim', [0 3.1], 'fontsize', 8);
cb = colorbar('peer', ax(end), 'location', 'east');
set(cb, 'position', [0.93 0.05 0.03 0.9/nens]);

