function avg = initarchivenc(vartable, Grd, Arch, filename)

% Start file

nc = netcdf(filename, 'clobber');  

% File description

nc.title = 'mixed_layer.m simulation output';
nc.history = sprintf('%s: %s', datestr(now), 'Mixed_layer simulation started');

% Add dimensions

nz = length(Grd.z);
nzp = length(Grd.zp);
nt = size(Arch.startdate, 1) + 1;

nc('depth') = nz;
nc('depth_edge') = nzp;
nc('time') = nt;

% Coordinate variables

nc{'depth'} = 'depth';
nc{'depth'}.unit = 'm';
nc{'depth'}(:) = Grd.z;

nc{'depth_edge'} = 'depth_edge';
nc{'depth_edge'}.unit = 'm';
nc{'depth_edge'}(:) = Grd.zp;

nc{'startdate'} = 'time';
nc{'startdate'}.unit = 'serial date number';
nc{'startdate'}(:) = datenum([Arch.startdate(1,:); Arch.startdate]);

nc{'enddate'} = 'time';
nc{'enddate'}.unit = 'serial date number';
nc{'enddate'}(:) = datenum([Arch.enddate(1,:); Arch.enddate]);

nc{'middate'} = 'time';
nc{'middate'}.unit = 'serial date number';
nc{'middate'}(:) = datenum([Arch.middate(1,:); Arch.middate]);

% Data variables

for ivar = 1:size(vartable,1)
    
    fprintf('  Initializing variable %d of %d\n', ivar, size(vartable,1));
    
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

% Set up average arrays

sz1 = cellfun(@size, vartable(:,1), 'uni', 0);
avg = cellfun(@zeros, sz1, 'uni', 0);


