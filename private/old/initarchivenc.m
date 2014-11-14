function avg = initarchivenc(vartable, Grd, Arch, filename)

% Start file

mode = bitor(nc_clobber_mode, nc_64bit_offset_mode);
nc_create_empty(filename, mode);

% File description

nc_attput(filename, nc_global, 'title', 'mixed_layer.m simulation output' );
nc_addhist(filename, 'Mixed_layer simulation started');

% Add dimensions

nz = length(Grd.z);
nzp = length(Grd.zp);
nt = size(Arch.startdate, 1) + 1;

nc_add_dimension(filename, 'depth', nz);
nc_add_dimension(filename, 'depth_edge', nzp);
nc_add_dimension(filename, 'time', nt);

% Coordinate variables

Var.Name = 'depth';
Var.Dimension = {'depth'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, Grd.z);
nc_attput(filename, Var.Name, 'unit', 'm');

Var.Name = 'depth_edge';
Var.Dimension = {'depth_edge'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, Grd.zp);
nc_attput(filename, Var.Name, 'unit', 'm');

Var.Name = 'startdate';
Var.Dimension = {'time'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, datenum([Arch.startdate(1,:); Arch.startdate]));
nc_attput(filename, Var.Name, 'unit', 'serial date number');

Var.Name = 'enddate';
Var.Dimension = {'time'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, datenum([Arch.enddate(1,:); Arch.enddate]));
nc_attput(filename, Var.Name, 'unit', 'serial date number');

Var.Name = 'middate';
Var.Dimension = {'time'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, datenum([Arch.middate(1,:); Arch.middate]));
nc_attput(filename, Var.Name, 'unit', 'serial date number');

% Data variables

for ivar = 1:size(vartable,1)
    
    Var.Name = vartable{ivar,2};
    
    if isequal(size(vartable{ivar,1}), [nz 1])
        Var.Dimension = {'time', 'depth'};
        start = [0 0];
        count = [1 nz];
    elseif isequal(size(vartable{ivar,1}), [nzp 1])
        Var.Dimension = {'time', 'depth_edge'};
        start = [0 0];
        count = [1 nzp];
    elseif isequal(size(vartable{ivar,1}), [1 1])
        Var.Dimension = {'time'};
        start = 0;
        count = 1;
    else
        error('Archive variable has unexpected dimensions');
    end
    
    nc_addvar(filename, Var);
    nc_varput(filename, Var.Name, vartable{ivar,1}, start, count);
    
    nc_attput(filename, Var.Name, 'long_name', vartable{ivar,3});
    nc_attput(filename, Var.Name, 'unit', vartable{ivar,4});
    
end

% Set up average arrays

sz1 = cellfun(@size, vartable(:,1), 'uni', 0);
avg = cellfun(@zeros, sz1, 'uni', 0);


