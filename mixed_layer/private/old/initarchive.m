function [avg, fid] = initarchive(data, ntime)

% Calculate size of data arrays

sz1 = cellfun(@size, data, 'uni', 0);
sz2 = cellfun(@(x) [x ntime], sz1, 'uni', 0);

% Set up average values arrays

avg = cellfun(@zeros, sz1, 'uni', 0);

% Open archive files

for ifile = 1:length(avg)
    filename = sprintf('mltemp%d.bin', ifile);
    fid(ifile) = fopen(filename, 'w');
end

% Record initial conditions (this is the only archive value that reflects
% only a single timestep rather than an averaged interval)

for idata = 1:length(data)
    fwrite(fid(idata), data{idata}, 'double');
end
