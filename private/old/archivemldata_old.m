 function avg = archivemldata(Grd, Arch, outfile, it, data)
%ARCHIVEMLDATA Write mixed_layer results to file
%
% avg = archivemldata(Grd, Arch, outfile, it, data)
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
%               or scalars.  Column 2-4 are strings with the variables'
%               short names, long names, and units, respectively.
%
% Output variables:
%
%   avg:        nvar x 1 cell array, holding the newly-averaged values each
%               variable.  Averages are performed over each archiving
%               period.

% Copyright 2010 Kelly Kearney



% Note: I use a mix of the netcdf toolbox and snctools because the netcdf
% method seems to be much faster to create the initial file, but much
% slower when writing data to the file as the simulation runs.

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

% Create new output file, if beginning of new file timestep

if Arch.isnewfile(it)
    
    % Create new file
    
    timeidx = Arch.filedates(Arch.fileidx(it,1),:);
    timeidx = timeidx(1):timeidx(2);
    sdate = Arch.startdate(timeidx,:); %TODO add initial time step?
    edate = Arch.enddate(timeidx,:);
    mdate = Arch.middate(timeidx,:);
    
    createfile(outfile, Arch.fileidx(it), sdate, edate, mdate, data, Grd.z, Grd.zp, timeidx);
    
end

% Save to file, if end of averaging time step
    
if Arch.islast(it)
    
    filename = regexprep(outfile, '\.nc', sprintf('_%03d.tmp.nc', Arch.fileidx(it)));
    bin = Arch.fileidx(it,2);
    
    % Write to file (netcdf toolbox method)
    
%     nc = netcdf(filename, 'write');  
%     for ivar = 1:size(data,1)
%         nc{data{ivar,2}}(bin,:) = avg{ivar}; % TODO add bin+1 if initial time step
%     end
        
    % Write to file (snctools method, much faster on my computer)
    
    for ivar = 1:size(data,1)
        if isscalar(avg{ivar})
            start = bin-1;
            count = 1;
        elseif isvector(avg{ivar});
            start = [bin-1 0];
            count = [1 length(avg{ivar})];
        end
        nc_varput(filename, data{ivar,2}, avg{ivar}, start, count);
    end
    
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    
end

% Combine files if reached last time step

if it == Grd.nt
    
    nfile = size(Arch.filedates,1);
    
    if nfile > 1
        cmd = sprintf('ncrcat -O %s_[0-9][0-9][0-9].tmp.nc %s', outfile(1:end-3), outfile);
        [s,r] = system(cmd);
        if s
            error('Error concatenating files:%s', r);
        end
        for ifl = 1:nfile
            delete(sprintf('%s_%03d.tmp.nc', outfile(1:end-3), ifl));
        end
    else
        movefile(sprintf('%s_001.tmp.nc', outfile(1:end-3)), outfile);
    end
    
end
    


function createfile(outfile, ifl, sdate, edate, mdate, vartable, z, zp, times)
%
% Old version, for use w/ Matlab 2007b or earlier w/ netcdf toolbox

filename = regexprep(outfile, '\.nc', sprintf('_%03d.tmp.nc', ifl));
% filename = sprintf('%s_%03d.tmp.nc', outfile, ifl);

% Start file

nc = netcdf(filename, 'clobber');  

% File description

% nc.title = 'mixed_layer.m simulation output';
% nc.history = sprintf('%s: %s', datestr(now), 'Mixed_layer simulation started');

% Add dimensions

nz = length(z);
nzp = length(zp);
nt = size(sdate,1);

nc('depth') = nz;
nc('depth_edge') = nzp;
nc('time') = 0;

% Coordinate variables

nc{'depth'} = 'depth';
nc{'depth'}.unit = 'm';
nc{'depth'}(:) = z;

nc{'depth_edge'} = 'depth_edge';
nc{'depth_edge'}.unit = 'm';
nc{'depth_edge'}(:) = zp;

nc{'startdate'} = 'time';
nc{'startdate'}.unit = 'serial date number';
nc{'startdate'}(1:nt) = datenum(sdate);

nc{'enddate'} = 'time';
nc{'enddate'}.unit = 'serial date number';
nc{'enddate'}(1:nt) = datenum(edate);

nc{'middate'} = 'time';
nc{'middate'}.unit = 'serial date number';
nc{'middate'}(1:nt) = datenum(mdate);

% nc{'time'} = 'time';
% nc{'time'}(:) = times;

% Data variables

for ivar = 1:size(vartable,1)
    
    if isequal(size(vartable{ivar,1}), [nz 1])
        nc{vartable{ivar,2}} = {'time', 'depth'};
    elseif isequal(size(vartable{ivar,1}), [nzp 1])
        nc{vartable{ivar,2}} = {'time', 'depth_edge'};
    elseif isequal(size(vartable{ivar,1}), [1 1])
        nc{vartable{ivar,2}} = 'time';
    else
        error('Archive variable has unexpected dimensions');
    end
    
    nc{vartable{ivar,2}}(1,:) = vartable{ivar,1};
    
    nc{vartable{ivar,2}}.long_name = vartable{ivar,3};
    nc{vartable{ivar,2}}.unit      = vartable{ivar,4};
    
end
close(nc);


function createfile(outfile, ifl, sdate, edate, mdate, vartable, z, zp, times)

filename = regexprep(outfile, '\.nc', sprintf('_%03d.tmp.nc', ifl));
% filename = sprintf('%s_%03d.tmp.nc', outfile, ifl);

% Start file

ncid = netcdf.create(filename, 'clobber');

% nc = netcdf(filename, 'clobber');  

% File description

% nc.title = 'mixed_layer.m simulation output';
% nc.history = sprintf('%s: %s', datestr(now), 'Mixed_layer simulation started');

% Add dimensions

nz = length(z);
nzp = length(zp);
nt = size(sdate,1);

nc('depth') = nz;
nc('depth_edge') = nzp;
nc('time') = 0;

% Coordinate variables

nc{'depth'} = 'depth';
nc{'depth'}.unit = 'm';
nc{'depth'}(:) = z;

nc{'depth_edge'} = 'depth_edge';
nc{'depth_edge'}.unit = 'm';
nc{'depth_edge'}(:) = zp;

nc{'startdate'} = 'time';
nc{'startdate'}.unit = 'serial date number';
nc{'startdate'}(1:nt) = datenum(sdate);

nc{'enddate'} = 'time';
nc{'enddate'}.unit = 'serial date number';
nc{'enddate'}(1:nt) = datenum(edate);

nc{'middate'} = 'time';
nc{'middate'}.unit = 'serial date number';
nc{'middate'}(1:nt) = datenum(mdate);

% nc{'time'} = 'time';
% nc{'time'}(:) = times;

% Data variables

for ivar = 1:size(vartable,1)
    
    if isequal(size(vartable{ivar,1}), [nz 1])
        nc{vartable{ivar,2}} = {'time', 'depth'};
    elseif isequal(size(vartable{ivar,1}), [nzp 1])
        nc{vartable{ivar,2}} = {'time', 'depth_edge'};
    elseif isequal(size(vartable{ivar,1}), [1 1])
        nc{vartable{ivar,2}} = 'time';
    else
        error('Archive variable has unexpected dimensions');
    end
    
    nc{vartable{ivar,2}}(1,:) = vartable{ivar,1};
    
    nc{vartable{ivar,2}}.long_name = vartable{ivar,3};
    nc{vartable{ivar,2}}.unit      = vartable{ivar,4};
    
end
close(nc);
