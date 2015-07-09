function postprocesswce(folder, filter, outfolder)
%POSTPROCESSWCE Reformat mixed_layer-wce output
%
% postprocesswce(folder, filter, outfolder)
%
% The mixed_layer model outputs one file per simulation, but when I run
% ensembles I often want to look at only the biological state variables,
% across all ensemble members.  This function extracts each variable into
% its own ensemble x depth x time file.  It also calculates critter x
% ensemble x time x depth biomass, net production, predatory mortality, and
% non-predatory mortality for the set of simulations (allbio.nc) and
% depth-integrated values of that (allbio_int.nc).
%
% Note: this does the same thing as postprocess.ncl, but that was getting
% very slow and hitting memory issues for some reason.
%
% Input variables:
%
%   folder:     name of folder where simulation output is located.  Files in
%               this folder are expected to be named sequentially.
%
%   filter:     dir-style filter to use on folder to acquire proper list of
%               files in a single ensemble. Use '*.nc' if no filter
%               required.
%
%   outfolder:  folder to place new .nc files.  A new folder with this name
%               will be created under the folder/postprocessed. It will be
%               created if it doesn't already exist.  If post-processed
%               files already exist in this folder, those variable files
%               will be skipped; delete old ones to rerun this script.
%               However, the bio-calc files are all or nothing; if either
%               is missing, both are recreated.

% Copyright 2014 Kelly Kearney


pattern = fullfile(folder, filter);
Files = dirfull(pattern);

if isempty(Files)
    error('No files found');
end

% Get list of variables

iscrashed = regexpfound({Files.name}, 'crashed');
idxfull = find(~iscrashed, 1);

Info = ncinfo(Files(idxfull).name);

[tf, loc] = ismember({'depth', 'depth_edge', 'time'}, {Info.Dimensions.Name});
if ~all(tf)
    error('Unexpected dimension names in file %s', Files(idxfull).name);
end

nt = Info.Dimensions(loc(3)).Length;
nz = Info.Dimensions(loc(1)).Length;
nens = length(Files);

newfolder = fullfile(folder, 'postprocessed', outfolder);
if ~exist(newfolder, 'dir')
    mkdir(newfolder);
end

% Categorize variables

dimvars =  {'depth', 'depth_edge', 'startdate', 'middate', 'enddate'};
[sname, lname] = readwcevarnames(Files(idxfull).name);

flx = regexp({Info.Variables.Name}, '(gpp|res|exx|gra|pre|exc|ege|mor|dec)_(\d+)_(\d+)', 'tokens', 'once');
vidx = find(~cellfun('isempty', flx));
flx = cat(1, flx{:});
fidx = cellfun(@str2num, flx(:,2:3));
flxname = {Info.Variables(vidx).Name};

bdiagvars = {...
'PSpsmax'   
'PLpsmax'   
'PSlightlim'
'PLlightlim'
'PSno3lim'  
'PLno3lim'  
'PSnh4lim'  
'PLnh4lim'  
'PLsilim'   
'PSfelim'   
'PLfelim'   
'I'         
'kappa'     
'kp'        
'PSfe2n'    
'PLfe2n'    
'PSfedef'   
'PLfedef'   
'extrasi'};


isdiag = ismember({Info.Variables.Name}, bdiagvars);

isbiovar = ismember({Info.Variables.Name}, [sname flxname]);
idx = find(isbiovar | isdiag);
nv = length(idx);

% Copy variable data into new files, one file per variable

fprintf('Collecting variables...\n');

cpb = ConsoleProgressBar();
cpb.setMinimum(0);      % minimum value of progress range [min max]
cpb.setMaximum(nens);   % maximum value of progress range [min max]
    

for ii = 1:nv
    
    iv = idx(ii);
    
    newfile = fullfile(newfolder, [Info.Variables(iv).Name '.nc']);
    if exist(newfile, 'file')
        continue
    end
    
    sz = [nens Info.Variables(iv).Size];

    % Create new file

    New.Name = '/';
    New.Dimensions = Info.Variables(iv).Dimensions;
    nd = length(New.Dimensions) + 1;
    New.Dimensions(nd).Name = 'ensemble';
    New.Dimensions(nd).Length = nens;
    New.Dimensions(nd).Unlimited = false;
    New.Dimensions = New.Dimensions([end 1:end-1]);

    New.Variables(1) = Info.Variables(iv);
    New.Variables(1).Dimensions = New.Dimensions;
    New.Variables(1).Size = sz;
    New.Variables(1).FillValue = -999;

    
    ncwriteschema(newfile, New);

    fprintf('  %d/%d: %s\n', ii, nv, Info.Variables(iv).Name);
    cpb.start();
    for ie = 1:nens
        cpb.setValue(ie);

        data = ncread(Files(ie).name, Info.Variables(iv).Name);
        ncwrite(newfile, Info.Variables(iv).Name, permute(data, [3 1 2]), [ie 1 1]);
    end
    cpb.stop();
    fprintf('\n');

end

% Combine biology data into single file


biofile = fullfile(newfolder, 'allbio.nc');
biointfile = fullfile(newfolder, 'allbio_int.nc');
runbio = ~exist(biofile) | ~exist(biointfile);

if runbio
    
    fprintf('Creating biology variables...\n');

    ng = length(sname);

    B.Name = '/';
    B.Dimensions = struct('Name',   {'critter', 'ensemble', 'time', 'depth'}, ...
                          'Length', {ng, nens, nt, nz}, ...
                          'Unlimited', {false, false, true, false});

    if exist(biofile, 'file')
        warning('File %s exists; overwriting', biofile);
        delete(biofile);
    end
    ncwriteschema(biofile, B);

    dims = {'critter', ng, 'ensemble', nens, 'time', nt, 'depth', nz};

    nccreate(biofile, 'bio',  'Dimensions', dims);
    nccreate(biofile, 'prod', 'Dimensions', dims);
    nccreate(biofile, 'pmor', 'Dimensions', dims);
    nccreate(biofile, 'omor', 'Dimensions', dims);

    ncwriteatt(biofile, 'bio',  'long_name', 'Biomass');
    ncwriteatt(biofile, 'bio',  'unit',      'mol/m^3');
    ncwriteatt(biofile, 'prod', 'long_name', 'Net production');
    ncwriteatt(biofile, 'prod', 'unit',      'mol/m^3/s');
    ncwriteatt(biofile, 'pmor', 'long_name', 'Predatory mortality');
    ncwriteatt(biofile, 'pmor', 'unit',      'mol/m^3/s');
    ncwriteatt(biofile, 'omor', 'long_name', 'Non-predatory mortality');
    ncwriteatt(biofile, 'omor', 'unit',      'mol/m^3/s');

    % And an integrated-biomass version


   
    zedge = ncread(Files(idxfull).name, 'depth_edge');
    dz = -diff(zedge);

    Bi.Name = '/';
    Bi.Dimensions = struct('Name',   {'critter', 'ensemble', 'time'}, ...
                          'Length', {ng, nens, nt}, ...
                          'Unlimited', {false, false, true});
    if exist(biointfile, 'file')
        warning('File %s exists; overwriting', biointfile);
        delete(biointfile);
    end
    ncwriteschema(biointfile, B);

    dims = {'critter', ng, 'ensemble', nens, 'time', nt};

    nccreate(biointfile, 'bio',  'Dimensions', dims);
    nccreate(biointfile, 'prod', 'Dimensions', dims);
    nccreate(biointfile, 'pmor', 'Dimensions', dims);
    nccreate(biointfile, 'omor', 'Dimensions', dims);

    ncwriteatt(biointfile, 'bio',  'long_name', 'Biomass');
    ncwriteatt(biointfile, 'bio',  'unit',      'mol/m^2');
    ncwriteatt(biointfile, 'prod', 'long_name', 'Net production');
    ncwriteatt(biointfile, 'prod', 'unit',      'mol/m^2/s');
    ncwriteatt(biointfile, 'pmor', 'long_name', 'Predatory mortality');
    ncwriteatt(biointfile, 'pmor', 'unit',      'mol/m^2/s');
    ncwriteatt(biointfile, 'omor', 'long_name', 'Non-predatory mortality');
    ncwriteatt(biointfile, 'omor', 'unit',      'mol/m^2/s');


    % Read biological data

    cpb.setMinimum(0);      % minimum value of progress range [min max]
    cpb.setMaximum(ng);   % maximum value of progress range [min max]

    cpb.start();
    for ig = 1:ng

        cpb.setValue(ig);

        % Biomass

        data = ncread(fullfile(newfolder, [sname{ig} '.nc']), sname{ig});
        datai = sum(bsxfun(@times, data, permute(dz, [2 1 3])), 2);

        ncwrite(biofile,    'bio', permute(data,  [4 1 3 2]), [ig 1 1 1]);
        ncwrite(biointfile, 'bio', permute(datai, [4 1 3 2]), [ig 1 1]);

        % Production

        issrc = ismember(flx(:,1), {'pre','gra','gpp'})          & fidx(:,2) == ig;
        issnk = ismember(flx(:,1), {'exx', 'res', 'exc', 'ege'}) & fidx(:,1) == ig;

        data = zeros(nens,nz,nt);
        for ii = find(issrc)'
            data = data + ncread(fullfile(newfolder, [flxname{ii} '.nc']), flxname{ii});
        end
        for ii = find(issnk)'
            data = data - ncread(fullfile(newfolder, [flxname{ii} '.nc']), flxname{ii});
        end
        datai = sum(bsxfun(@times, data, permute(dz, [2 1 3])), 2);

        ncwrite(biofile,    'prod', permute(data,  [4 1 3 2]), [ig 1 1 1]);
        ncwrite(biointfile, 'prod', permute(datai, [4 1 3 2]), [ig 1 1]);

        % Predatory mortality

        ispmr = ismember(flx(:,1), {'pre','gra'})                & fidx(:,1) == ig;

        data = zeros(nens,nz,nt);
        for ii = find(ispmr)'
            data = data + ncread(fullfile(newfolder, [flxname{ii} '.nc']), flxname{ii});
        end
        datai = sum(bsxfun(@times, data, permute(dz, [2 1 3])), 2);

        ncwrite(biofile,    'pmor', permute(data,  [4 1 3 2]), [ig 1 1 1]);
        ncwrite(biointfile, 'pmor', permute(datai, [4 1 3 2]), [ig 1 1]);

        % Non-predatory mortality

        isomr = ismember(flx(:,1), {'mor'})                      & fidx(:,1) == ig;

        data = zeros(nens,nz,nt);
        for ii = find(isomr)'
            data = data + ncread(fullfile(newfolder, [flxname{ii} '.nc']), flxname{ii});
        end
        datai = sum(bsxfun(@times, data, permute(dz, [2 1 3])), 2);

        ncwrite(biofile,    'omor', permute(data,  [4 1 3 2]), [ig 1 1 1]);
        ncwrite(biointfile, 'omor', permute(datai, [4 1 3 2]), [ig 1 1]);
    end
    cpb.stop();
    fprintf('\n');
end


                      
















