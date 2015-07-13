function avg = archivedatanc(frac, islast, bin, avg, filename, data)

ndata = length(data);

% Do the averaging

for idata = 1:ndata
    avg{idata} = avg{idata} + frac .* data{idata};
end

% Archiving

if islast
   
    % Write to file
    
    for ivar = 1:size(data,1)
        if isscalar(avg{ivar})
            start = bin;
            count = 1;
        elseif isvector(avg{ivar});
            start = [bin 0];
            count = [1 length(avg{ivar})];
        end
        nc_varput(filename, data{ivar,2}, avg{ivar}, start, count);
    end
        
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end