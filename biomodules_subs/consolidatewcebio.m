function consolidatewcebio(folder)
%CONSOLIDATEWCEBIO Consolidate biological variables for wce run
%
% consolidatewcebio(folder)
%
% This function creates two new netcdf files that hold biomass, net
% production, predatory mortality, and non-predatory mortality for all
% biological state variables in a wce run.  Both files are added to the
% simulation folder.  The allbio.nc file includes depth x time x ensemble x
% critter (dimensions in reverse order when read via ncread) arrays, while
% allbio_int.nc holds depth-integrated time x ensemble x critter arrays.
% See critternames.mat for the order of the critter dimension.

% Copyright 2014 Kelly Kearney

if ~exist(fullfile(folder, 'temp.nc'), 'file')
    error('Postprocessed files not found');
end

Tinfo = ncinfo(fullfile(folder, 'temp.nc'));
        
[tf, loc] = ismember({'ensemble', 'time', 'depth'}, {Tinfo.Dimensions.Name});
nens = Tinfo.Dimensions(loc(1)).Length;
nt   = Tinfo.Dimensions(loc(2)).Length;
nz   = Tinfo.Dimensions(loc(3)).Length;


V = dir(fullfile(folder, '*.nc'));
isvar = ~strcmp({V.name}, 'dimensions.nc') & ~regexpfound({V.name}, 'sim\d\d\d\d.nc');
vars = regexprep({V(isvar).name}, '\.nc', '');
if isempty(vars)
    warning('Cannot consolidate biology; files not found');
    return
end

biofile = fullfile(folder, 'allbio.nc');
biointfile = fullfile(folder, 'allbio_int.nc');

critterfile = fullfile(folder, 'crittername.mat');
samplefile  = fullfile(folder, 'sim0001.nc');
dimfile     = fullfile(folder, 'dimensions.nc');

% Easier to just save the critter names to a .mat file than
% deal with netcdf character array.  For record, though, add
% names as an attribute to the dimensions file.

if exist(critterfile, 'file')
    N = load(critterfile);
    sname = N.sname; 
else
    error('Cannot consolidate biology; need critter names');
end

% Parse flux variables

flx = regexp(vars, '(gpp|res|exx|gra|pre|exc|ege|mor|dec)_(\d+)_(\d+)', 'tokens', 'once');
vidx = find(~cellfun('isempty', flx));
flx = cat(1, flx{:});
fidx = cellfun(@str2num, flx(:,2:3));
flxname = vars(vidx);

% Create new biology file and integrated-biology file

biofile = fullfile(folder, 'allbio.nc');
biointfile = fullfile(folder, 'allbio_int.nc');

ng = length(sname);

Dim(1).Name = 'critter';
Dim(1).Length = ng;
Dim(1).Unlimited = false;

Dim(2).Name = 'ensemble';
Dim(2).Length = nens;
Dim(2).Unlimited = false;

Dim(3).Name = 'time';
Dim(3).Length = nt;
Dim(3).Unlimited = true;

Dim(4).Name = 'depth';
Dim(4).Length = nz;
Dim(4).Unlimited = false;

B = struct('Name', '/', 'Format', 'netcdf4');

B.Variables(1).Name = 'bio';
B.Variables(1).Dimensions = Dim;
B.Variables(1).Attributes = struct(...
    'Name',  {'long_name', 'unit'}, ...
    'Value', {'Biomass',   'mol/m^3'});
B.Variables(1).Datatype = 'double';
B.Variables(1).FillValue = -9999;

B.Variables(2).Name = 'prod';
B.Variables(2).Dimensions = Dim;
B.Variables(2).Attributes = struct(...
    'Name',  {'long_name',      'unit'}, ...
    'Value', {'Net production', 'mol/m^3'});
B.Variables(2).Datatype = 'double';
B.Variables(2).FillValue = -9999;

B.Variables(3).Name = 'pmor';
B.Variables(3).Dimensions = Dim;
B.Variables(3).Attributes = struct(...
    'Name',  {'long_name',           'unit'}, ...
    'Value', {'Predatory mortality', 'mol/m^3'});
B.Variables(3).Datatype = 'double';
B.Variables(3).FillValue = -9999;

B.Variables(4).Name = 'omor';
B.Variables(4).Dimensions = Dim;
B.Variables(4).Attributes = struct(...
    'Name',  {'long_name',               'unit'}, ...
    'Value', {'Non-predatory mortality', 'mol/m^3'});
B.Variables(4).Datatype = 'double';
B.Variables(4).FillValue = -9999;

if exist(biofile, 'file')
    warning('File %s exists; overwriting', biofile);
    delete(biofile);
end
ncwriteschema(biofile, B);

Bi = B;
[Bi.Variables.Dimensions] = deal(Dim(1:3));
for iv = 1:length(Bi.Variables)
    Bi.Variables(iv).Attributes(2).Value = 'mol/m^2';
end

if exist(biointfile, 'file')
    warning('File %s exists; overwriting', biointfile);
    delete(biointfile);
end
ncwriteschema(biointfile, Bi);

zedge = ncread(dimfile, 'depth_edge');
dz = -diff(zedge);

% Read and rewrite biology data

fprintf('Collecting biology: %s\n', folder);
cpb = ConsoleProgressBar();
cpb.setMinimum(0);      
cpb.setMaximum(ng);   

cpb.start();
for ig = 1:ng

    cpb.setValue(ig);
    cpb.setText(sprintf('%2d/%d', ig, ng));

    % Biomass

    data = ncread(fullfile(folder, [sname{ig} '.nc']), sname{ig});
    datai = sum(bsxfun(@times, data, dz), 1);

    ncwrite(biofile,    'bio', permute(data,  [4 3 2 1]), [ig 1 1 1]);
    ncwrite(biointfile, 'bio', permute(datai, [4 3 2 1]), [ig 1 1]);

    % Production

    issrc = ismember(flx(:,1), {'pre','gra','gpp'})          & fidx(:,2) == ig;
    issnk = ismember(flx(:,1), {'exx', 'res', 'exc', 'ege'}) & fidx(:,1) == ig;

    data = zeros(nz,nt,nens);
    for ii = find(issrc)'
        data = data + ncread(fullfile(folder, [flxname{ii} '.nc']), flxname{ii});
    end
    for ii = find(issnk)'
        data = data - ncread(fullfile(folder, [flxname{ii} '.nc']), flxname{ii});
    end
    datai = sum(bsxfun(@times, data, dz), 1);

    ncwrite(biofile,    'prod', permute(data,  [4 3 2 1]), [ig 1 1 1]);
    ncwrite(biointfile, 'prod', permute(datai, [4 3 2 1]), [ig 1 1]);

    % Predatory mortality

    ispmr = ismember(flx(:,1), {'pre','gra'})                & fidx(:,1) == ig;

    data = zeros(nz,nt,nens);
    for ii = find(ispmr)'
        data = data + ncread(fullfile(folder, [flxname{ii} '.nc']), flxname{ii});
    end
    datai = sum(bsxfun(@times, data, dz), 1);

    ncwrite(biofile,    'pmor', permute(data,  [4 3 2 1]), [ig 1 1 1]);
    ncwrite(biointfile, 'pmor', permute(datai, [4 3 2 1]), [ig 1 1]);

    % Non-predatory mortality

    isomr = ismember(flx(:,1), {'mor'})                      & fidx(:,1) == ig;

    data = zeros(nz,nt,nens);
    for ii = find(isomr)'
        data = data + ncread(fullfile(folder, [flxname{ii} '.nc']), flxname{ii});
    end
    datai = sum(bsxfun(@times, data, dz), 1);

    ncwrite(biofile,    'omor', permute(data,  [4 3 2 1]), [ig 1 1 1]);
    ncwrite(biointfile, 'omor', permute(datai, [4 3 2 1]), [ig 1 1]);
end
cpb.stop();
fprintf('\n');