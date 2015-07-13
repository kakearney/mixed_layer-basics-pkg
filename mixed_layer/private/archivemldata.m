 function [avg, ncid, vid] = archivemldata(Grd, Arch, it, data, it0)
%ARCHIVEMLDATA Write mixed_layer results to file
%
% avg = archivemldata(Grd, Arch, outfile, it, data)
% [avg, fid] = archivemldata(Grd, Arch, outfile, it, data)
%
% New version should be compatible with all versions of Matlab, as long as
% snctools is installed and mexnc points to the proper backend (either
% newer internal stuff or the proper mex version).
%
% Input variables:
%
%   Grd:        1 x 1 structure holding spatial and temporal grid data for
%               mixed_layer simulation 
%
%   Arch:       1 x 1 structure holding data related to archiving
%
%   outfile:    name of output file
%
%   it:         index of current model time step
%
%   data:       nvar x 4 cell array of data to be written to file.  Column
%               1 holds the actual data values for the current model time
%               step, which must be either nz x 1 vectors, nzp x 1 vectors,
%               or scalars.  Columns 2-4 are strings with the variables'
%               short names, long names, and units, respectively.
%
%   it0:        index of first time step being simulated
%
% Output variables:
%
%   avg:        nvar x 1 cell array, holding the newly-averaged values each
%               variable.  Averages are performed over each archiving
%               period.

% Copyright 2010 Kelly Kearney


% If very first step, set up average arrays

if it == 1
    sz1 = cellfun(@size, data(:,1), 'uni', 0);
    avg = cellfun(@zeros, sz1, 'uni', 0);
else
    avg = Arch.avg;
end

% Include the current time step in the averages

ndata = size(data,1);
for idata = 1:ndata
    avg{idata} = avg{idata} + Arch.fraction(it) .* data{idata};
end

% Create output files

if it == it0
    
        [ncid, vid] = createoutputfiles(...
            Arch.outfile, ...
            Grd.z, ...
            Grd.zp, ...
            Arch.dateedge(1:end-1), ...
            Arch.dateedge(2:end), ...
            Arch.middate, ...
            data, ...
            Arch.iens);
else
    ncid = Arch.ncid;
    vid = Arch.vid;
end

% Write data to files

if Arch.islast(it)
    
    isd = cellfun(@(x) isequal(size(Grd.z), size(x)), avg);
    isde = cellfun(@(x) isequal(size(Grd.zp), size(x)), avg);
    iss = cellfun(@(x) isscalar(x), avg);
    
    start = cell(size(avg));
    [start{isd | isde}] = deal([0 Arch.bin(it)-1]);
    [start{iss}] = deal([Arch.bin(it)-1]);
    
    count = cell(size(avg));
    [count{isd}] = deal([length(Grd.z) 1]);
    [count{isde}] = deal([length(Grd.zp) 1]);
    [count{iss}] = deal([1]);
    
    for iv = 1:length(avg)
        netcdf.putVar(ncid, vid(iv), start{iv}, count{iv}, avg{iv});
    end
    
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end


