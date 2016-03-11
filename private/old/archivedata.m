function avg = archivedata(frac, islast, bin, avg, fid, data)

ndata = length(data);

% Do the averaging

for idata = 1:ndata
    avg{idata} = avg{idata} + frac .* data{idata};
end

% Archiving

if islast
   
    % Transfer to output
    
    for idata = 1:ndata
        fwrite(fid(idata), avg{idata}, 'double');
    end
    
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end