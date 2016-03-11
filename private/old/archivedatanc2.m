function avg = archivedatanc(frac, islast, bin, avg, filename, data)

ndata = length(data);

% Do the averaging

for idata = 1:ndata
    avg{idata} = avg{idata} + frac .* data{idata};
end

% Archiving

if islast
    
    nc = netcdf(filename, 'write');  
   
    % Write to file
    
    for ivar = 1:size(data,1)
        nc{data{ivar,2}}(bin+1,:) = avg{ivar};
    end
        
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end