function recovercrashed(varargin)
%RECOVERCRASHED Create netcdf output for a partially-run mixed_layer sim
% 
% recovercrashed(varargin)
%
% This function recovers the data from a crashed mixed_layer run.  It
% assumes that the crashed one was the most recently run mixed_layer
% simulation, and that the mltemp.bin and mltemp.mat files in the current
% directory are associated with it.
%
% Run with the same input as was used for mixed_layer.

% Copyright 2010 Kelly Kearney

% Run initilization-only to get all run paramters

P = mixed_layer(varargin{:}, 'stopafterinit', true);

nfile = length(P.Arch);

% Figure out which temp file corresponds to the run

AllFiles = dir(fullfile(P.In.tempdir, '*.mat'));
vars = arrayfun(@(X) who('-file', fullfile(P.In.tempdir,X.name)), AllFiles, 'uni', 0);
ismlmat = cellfun(@(x) any(strcmp(x, 'mlmarker')), vars);
AllFiles = AllFiles(ismlmat);
Tmp = arrayfun(@(X) load(fullfile(P.In.tempdir,X.name), 'outfile'), AllFiles); 
filematch = {Tmp.outfile}';

for ifl = 1:nfile
    
    ismatch = strcmp(filematch, P.Arch(ifl).outfile);
    
%     ismatch = false(size(Files));
%     for ii = 1:length(Files)
%         cfile = fullfile(tempdir, Files(ii).name);
%         vars = who('-file', cfile);
%         if ismember('mlmarker', vars)
%             Tmp = load(cfile, 'mlmarker', 'outfile');
%             ismatch(ii) = strcmp(Tmp.mlmarker, 'mixed_layer temp file') && strcmp(Tmp.outfile, P.Arch(ifl).outfile);
%         end
%     end

    Files = AllFiles(ismatch);
    nc = length(Files);

    if nc < 1
        error('MIXED_LAYER:norecoverfile', 'No mixed_layer temp files located that match this input');
    elseif nc > 1
        fprintf('Multiple temp files found that match this input:\n');
        for ic = 1:nc
            fprintf('  %d: %s (%s)\n', ic, Files(ic).name, Files(ic).date);
        end
        idx = input('Which file do you want to use? ');
        matfile = fullfile(P.In.tempdir, Files(idx).name);
    else
        matfile = fullfile(P.In.tempdir, Files.name);
    end
    binfile = regexprep(matfile, '.mat', '.bin');

    Info = dir(binfile);
    if Info.bytes == 0
        continue
%         error('Binary file empty');
    end

    A = load(matfile);

    % Convert temporary file to netcdf

    outfile = regexprep(P.Arch(ifl).outfile, '\.nc', '.crashed.nc');

    sdate = P.Arch(ifl).dateedge(1:end-1);
    edate = P.Arch(ifl).dateedge(2:end);

    binary2netcdf(binfile, outfile, sdate, edate, ...
                  P.Arch(ifl).middate, P.Grd.z, P.Grd.zp, A.data, ...
                  P.Arch(ifl).nper, P.In.verbose);
end
