function binary2netcdf(binfile, ncfile, sdate, edate, mdate, z, zp, vartable, ntperread)

sz = cellfun(@(x) size(x,1), vartable(:,1));
npertime = sum(sz);
nvar = length(sz);

% How many time steps were written to this file?

fid = fopen('mltemp.bin');
fseek(fid, 0, 'eof');
bits = ftell(fid);

nt = bits/(8.*npertime);
if nt < size(sdate,1)
    warning('Binary file does not include full simulation');
    sdate = sdate(1:nt,:);
    edate = edate(1:nt,:);
    mdate = mdate(1:nt,:);
end

frewind(fid);

%--------------------------
% Create netcdf file
%--------------------------

% Start file

mode = bitor(nc_clobber_mode, nc_64bit_offset_mode);
nc_create_empty(ncfile, mode);


% Add dimensions

nz = length(z);
nzp = length(zp);

nc_add_dimension(ncfile, 'depth', nz);
nc_add_dimension(ncfile, 'depth_edge', nzp);
nc_add_dimension(ncfile, 'time', 0);

% Coordinate variables

Var.Name = 'depth';
Var.Dimension = {'depth'};
nc_addvar(ncfile, Var);
nc_varput(ncfile, Var.Name, z);
nc_attput(ncfile, Var.Name, 'unit', 'm');

Var.Name = 'depth_edge';
Var.Dimension = {'depth_edge'};
nc_addvar(ncfile, Var);
nc_varput(ncfile, Var.Name, zp);
nc_attput(ncfile, Var.Name, 'unit', 'm');

Var.Name = 'startdate';
Var.Dimension = {'time'};
nc_addvar(ncfile, Var);
nc_varput(ncfile, Var.Name, datenum(sdate));
nc_attput(ncfile, Var.Name, 'unit', 'serial date number');

Var.Name = 'enddate';
Var.Dimension = {'time'};
nc_addvar(ncfile, Var);
nc_varput(ncfile, Var.Name, datenum(edate));
nc_attput(ncfile, Var.Name, 'unit', 'serial date number');

Var.Name = 'middate';
Var.Dimension = {'time'};
nc_addvar(ncfile, Var);
nc_varput(ncfile, Var.Name, datenum(mdate));
nc_attput(ncfile, Var.Name, 'unit', 'serial date number');

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
    
    nc_addvar(ncfile, Var);
    
    nc_attput(ncfile, Var.Name, 'long_name', vartable{ivar,3});
    nc_attput(ncfile, Var.Name, 'unit', vartable{ivar,4});
    
end


%--------------------------
% Transfer data from binary
% file to netcdf file
%--------------------------

tstart = 0:ntperread:(nt-1);
tcount = ones(size(tstart)).*ntperread;

if mod(nt, ntperread)
    tcount(end) = mod(nt, ntperread);
end


for ii = 1:length(tstart)

    alldata = fread(fid, npertime.*tcount(ii), 'double');
    alldata = reshape(alldata, npertime, tcount(ii));
    alldata = mat2cell(alldata, sz, tcount(ii));

    for iv = 1:nvar

        if sz(iv) == 1
            start = tstart(ii);
            count = tcount(ii);
        else
            start = [tstart(ii) 0];
            count = [tcount(ii) sz(iv)];
        end

        nc_varput(ncfile, vartable{iv,2}, alldata{iv}', start, count);

    end
   
    
end

