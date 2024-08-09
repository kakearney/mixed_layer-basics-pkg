function runmixedlayer(In, varargin)
%RUNMIXEDLAYER Run mixed_layer model for multiple inputs
%
% runmixedlayer(In, p1, v1, ...)
%
% This function runs mixed_layer for multiple input sets.  It provides a
% log file to note any errors.  It will skip over any sets where it finds
% an existing file of the specified name, so can be used to complete
% partly-run ensembles. 
%
% Input variables:
%
%   In:         n x 1 structure holding mixed_layer input fields.  Each
%               structure element corresponds to one set of input
%               variables.
%
% Optional input variables (passed as parameter/value pairs):
%
%   name:       base name for run.  A folder with this name (or this name
%               modified by outputextension of the runs) will be placed
%               under the folder specified below, and will contain the
%               logfile name.log as well as all .nc files produced by the
%               run. ['test']
%
%   folder:     location where new results folder will be placed (new
%               folder will will be named <folder>/<name_outputextension>
%               ['.']
%
%   usepar:     true to run simulations in parallel via parpool and
%               parfor, false to run one at a time. [false]
%
%   nlabs:      number of labs (i.e. workers) to run if running in parallel
%               [4] 
%
%   idx:        indices of input structure to actually run after setup
%               [1:length(In)]
%
%   idxmissing: indices of ensemble members expected to be missing when
%               postprocessing begins, to be replaced with all-fill-value
%               files (will only occur if model crashes before first
%               archiving step; should only happen if hot-start file is
%               missing).
%
%   postprocess: logical scalar indicating whether to postprocess the
%               output files, converting from one file per ensemble member
%               to one file per variable. [false]
%
%   consolidatebio: logical scalar indicating whether to run
%               consolidatwcebio (only applicable to wce runs). [false]
%
%   p#/v#:      any additional mixed_layer input variables that will be
%               applied to all runs, passed as parameter/value pairs

% Copyright 2012-2014 Kelly Kearney

%------------------------------
% Setup
%------------------------------

% Parse options

Opt.folder = '.';
Opt.usepar = false;
Opt.nlabs = 4;
Opt.idx = 1:numel(In);
Opt.name = 'test';
Opt.postprocess = false;
Opt.consolidatebio = false;
Opt.idxmissing = [];

[Opt, pv] = parsepv(Opt, varargin, 'returnextra');

nens = numel(In);
In = In(:);

Opt.idxmissing = reshape(Opt.idxmissing, 1, []);

% Create folder, if necessary

if ~exist(Opt.folder, 'dir')
    mkdir(Opt.folder);
end

% Check if any ensemble members have already been run

fields = {'tarch', 'beginarchive', 'endarchive'};
pvtmp = reshape(pv,2,[]);
Pv = cell2struct(pvtmp(2,:)', pvtmp(1,:)');

if isfield(Pv, 'verbose')
    verb = true(nens,1) .* Pv.verbose;
else
    verb = true(nens,1);
    if isfield(In, 'verbose')
        isemp = cellfun('isempty', {In.verbose});
        verb(~isemp) = [In(~isemp).verbose];
    end
end

S = warning('off', 'all'); % suppress duplicate fieldnames warning, untagged
Tmp = mergestruct(In(1), Pv); 
warning(S);

if isfield(Tmp, 'outputextension') && ~isempty(Tmp.outputextension)
    simfolder = cellfun(@(x) fullfile(Opt.folder,[Opt.name '_' x]), Tmp.outputextension, 'uni', 0);
else
    isfld = isfield(Tmp, fields);
    ntmp = ones(1,3);
    ntmp(isfld) = cellfun(@(x) length(Tmp.(x)), fields(isfld));
    ntmp = max(ntmp);
    if ntmp == 1
        simfolder = {fullfile(Opt.folder, Opt.name)};
    else
        simfolder = arrayfun(@(x) fullfile(Opt.folder,sprintf('%s_%d',Opt.name,x)), (1:ntmp)', 'uni', 0);
    end
end
    
Files = cellfun(@(x) dir(fullfile(x,'sim*.nc')), simfolder, 'uni', 0);
idx = cellfun(@(X) cellfun(@(x) str2double(x{1}), regexp({X.name}, 'sim(\d+)', 'tokens', 'once')), Files, 'uni', 0);

% If I find the final post-processed files (using temp.nc as the marker,
% since it's always there regardless of the bio module), I'm going to
% assume things have been completed and the intermediate files
% intentionally deleted.  If you want to rerun, delete or move these files.

for is = 1:length(simfolder)
    if exist(fullfile(simfolder{is}, 'temp.nc'))
        idx{is} = 1:nens;
    end
end

% Indices to run

idx = cellfun(@(x) x(:), idx, 'uni', 0);
idx = unique(cat(1,idx{:})); % Assuming here that any mismatch between folders is intentional

Opt.idx = setdiff(Opt.idx, idx);
runflag = ismember(1:length(In), Opt.idx); %b/c parfor needs consecutive

%------------------------------
% Run simulations
%------------------------------

% Log file

logfile = fullfile(Opt.folder, [Opt.name datestr(now, '_yyyymmddHHMM') '.log']);

% Set up parallel, if necessary

pcloseflag = false;
p = gcp('nocreate');
if Opt.usepar
    if isempty(p) && any(runflag)
        p = parpool(Opt.nlabs);
        pcloseflag = true;
    end
    if isempty(p)
        nworker = 0;
    else
        nworker = p.NumWorkers;
    end
else
    nworker = 0;
end
    
% Run simulations

parfor (iw = 1:nens, nworker)
    
    if runflag(iw)
    
        file = fullfile(Opt.folder, Opt.name);
        if ~verb(iw)
            fprintf('%s: %s: %d\n', datestr(now), file, iw);
        end

        try
            mixed_layer(file, In(iw), pv{:}, 'iens', iw);
        catch Me
            stackdata = [{Me.stack.file}; {Me.stack.line}];

            fid = fopen(logfile, 'at');

            fprintf(fid, '----------------------------\nDate:\t%s\nRun:\t%d\n\n', datestr(now), iw);
            fprintf(fid, '%s\n', Me.message);
            fprintf(fid, '  File: %s, line %d\n', stackdata{:});
            fprintf(fid, '\n');
            fclose(fid);


        end
    end
    
end

if pcloseflag
    delete(p);
end

%------------------------------
% Postprocess: 1 file per sim 
% with all vars to 1 file per
% var with all ensemble members
%------------------------------

if Opt.postprocess

    for is = 1:length(simfolder)
        
        fprintf('Postprocessing %s\n', simfolder{is});
        
        % Make sure all ensemble members have been run; makes no sense to
        % postprocess before this is done (if postprocessing has already
        % been completed, this will skip that as well)
        
        Tmp = dir(fullfile(simfolder{is}, 'sim*.nc'));
        for ii = 1:length(Tmp)
            Tmp(ii).name = fullfile(Tmp(ii).folder, Tmp(ii).name);
        end
        
        if isempty(Tmp) && exist(fullfile(simfolder{is}, 'temp.nc'), 'file')
            warning('Folder %s already post-processed', simfolder{is});
            continue % already post-processed (if started from crash)
        end
        
        dlen = zeros(length(Tmp),1);
        for itmp = 1:length(Tmp)
            ncid = netcdf.open(Tmp(itmp).name);
            dmid = netcdf.inqDimID(ncid, 'time');
            [blah, dlen(itmp)] = netcdf.inqDim(ncid, dmid);
            netcdf.close(ncid);
        end
        [maxlen, idxfull] = max(dlen);
        
        % Some files may be intentionally missing (for example, sims
        % hot-starting from an ensemble with a few crashed members).
        % Create empty files for these
        
        Info = ncinfo(Tmp(idxfull).name);
        Info = rmfield(Info, 'Filename');
        for im = Opt.idxmissing
            newfile = fullfile(simfolder{is}, sprintf('sim%04d.nc', im));
            if exist(newfile, 'file')
               delete(newfile); % If crashed, start over, just in case 
            end
            ncwriteschema(newfile, Info);
        end
        
        % Figure out if there are any files missing that weren't expected
        % to be
        
        Tmp = dirfull(fullfile(simfolder{is}, 'sim*.nc'));
         
        idx = cellfun(@(x) str2double(x{1}), regexp({Tmp.name}, 'sim(\d+)', 'tokens', 'once'));
        idxmiss = setdiff(1:nens, idx);
        if ~isempty(idxmiss)
            str = sprintf('%d,', idxmiss);
            str = str(1:end-1);
            warning('Missing output file in %s (%s); cannot post-process', simfolder{is}, str);
            continue
        end
        
        % Check if files ran to completion or not
        
        len = nan(nens,1);
        for ii = 1:nens
            ncid = netcdf.open(Tmp(ii).name);
            dmid = netcdf.inqDimID(ncid, 'time');
            [blah, len(ii)] = netcdf.inqDim(ncid, dmid);
            netcdf.close(ncid);
        end
        nt = max(len);
        ispartial = len < nt;

        Info = ncinfo(Tmp(1).name);
        vars = {Info.Variables.Name}; 
        
        % Increase time dimension where necessary

        idxpart = find(ispartial);
        for ip = idxpart'
            ncwrite(Tmp(ip).name, 'temp', NaN, [1 nt]);            
        end

        % Concatenate files for each variable
             
        fprintf('Concatenating ensembles\n');
        cpb = ConsoleProgressBar();
        cpb.setMinimum(0);      
        cpb.setMaximum(length(vars));
        
        cpb.start();
        
        for iv = 1:length(vars)
            
            cpb.setValue(iv);
            str = sprintf('%3d of %3d: %s', iv, length(vars), vars{iv});
            cpb.setText(str);
            
            outfile = fullfile(simfolder{is}, [vars{iv} '.nc']);
            if ~exist(outfile, 'file')
                cmd = sprintf('ncecat -v %s -u ensemble -n %d,4,1 %s %s', vars{iv}, nens, Tmp(1).name, outfile);
                [s,r] = system(cmd);
                if s
                    error(r);
                end
            end
        end 
        cpb.stop();
        fprintf('\n');
        
        % For wce runs, need to preserve critter name order
        
        if isfield(In, 'biofun') && strcmp(func2str(In(1).biofun), 'wce')
            critterfile = fullfile(simfolder{is}, 'crittername.mat');

            [sname, lname] = readwcevarnames(Tmp(1).name);
            save(critterfile, 'sname', 'lname');
    
            critterstr = sprintf('%s,', sname{:});
            dimfile = fullfile(simfolder{is}, 'dimensions.nc');
            ncwriteatt(dimfile, '/', 'critters', critterstr(1:end-1));
            
        end
        
        % At this point, we don't need the original files anymore; all data
        % has been replicated elsewhere.

        arrayfun(@(X) delete(X.name), Tmp);
    end
end

%------------------------------
% Consolidate biological 
% variables
%------------------------------

if Opt.consolidatebio && isfield(In, 'biofun') 
    
    switch func2str(In(1).biofun)
        case 'wce'
    
            for is = 1:length(simfolder)
                consolidatewcebio(simfolder{is});
            end
            
        case 'nemurokak'
            for is = 1:length(simfolder)
                consolidatenemurokakbio(simfolder{is});
            end
        otherwise
            warning('consolidatebio only an option for wce and nemurokak runs');
    end
           
    
end







