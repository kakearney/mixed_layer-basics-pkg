function binary2netcdf(binfile, ncfile, sdate, edate, mdate, z, zp, vartable, ntperread)
%BINARY2NETCDF
%
% binary2netcdf(binfile, ncfile, sdate, edate, mdate, z, zp, vartable, ...
%               ntperread)
%
% This function creates the netcdf output file for mixed_layer.  It
% requires Matlab R2009a or later, which has netcdf support.
% binary2netcdf_snc can replace this function for older Matlab versions,
% but that version is currently absurdly slow.

% Copyright 2010 Kelly Kearney

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

fprintf('  Creating file\n');

% Start file

mode = bitor(netcdf.getConstant('NC_CLOBBER'), netcdf.getConstant('NC_64BIT_OFFSET')); 
ncid = netcdf.create(ncfile, mode);

% Add dimensions

nz = length(z);
nzp = length(zp);

zid = netcdf.defDim(ncid, 'depth', nz);
zpid = netcdf.defDim(ncid, 'depth_edge', nzp);
tid = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));

% Coordinate variables

cvid(1) = netcdf.defVar(ncid, 'depth',      'double', zid);
cvid(2) = netcdf.defVar(ncid, 'depth_edge', 'double', zpid);
cvid(3) = netcdf.defVar(ncid, 'startdate',  'double', tid);
cvid(4) = netcdf.defVar(ncid, 'enddate',    'double', tid);
cvid(5) = netcdf.defVar(ncid, 'middate',    'double', tid);

netcdf.putAtt(ncid, cvid(1), 'unit', 'm');
netcdf.putAtt(ncid, cvid(2), 'unit', 'm');
netcdf.putAtt(ncid, cvid(3), 'unit', 'serial date number');
netcdf.putAtt(ncid, cvid(4), 'unit', 'serial date number');
netcdf.putAtt(ncid, cvid(5), 'unit', 'serial date number');

netcdf.endDef(ncid);

netcdf.putVar(ncid, cvid(1), z);
netcdf.putVar(ncid, cvid(2), zp);
netcdf.putVar(ncid, cvid(3), 0, nt, datenum(sdate));
netcdf.putVar(ncid, cvid(4), datenum(edate));
netcdf.putVar(ncid, cvid(5), datenum(mdate));

% Data variables

netcdf.reDef(ncid);

for ivar = 1:size(vartable,1)
    
    if isequal(size(vartable{ivar,1}), [nz 1])
        vid(ivar) = netcdf.defVar(ncid, vartable{ivar,2}, 'double', [zid tid]);
    elseif isequal(size(vartable{ivar,1}), [nzp 1])
        vid(ivar) = netcdf.defVar(ncid, vartable{ivar,2}, 'double', [zpid tid]);
    elseif isequal(size(vartable{ivar,1}), [1 1])
        vid(ivar) = netcdf.defVar(ncid, vartable{ivar,2}, 'double', tid);
    end
    
    netcdf.putAtt(ncid, vid(ivar), 'long_name', vartable{ivar,3});
    netcdf.putAtt(ncid, vid(ivar), 'unit',      vartable{ivar,4});
    
end

netcdf.endDef(ncid);

%--------------------------
% Transfer data from binary
% file to netcdf file
%--------------------------

tstart = 0:ntperread:(nt-1);
tcount = ones(size(tstart)).*ntperread;

if mod(nt, ntperread)
    tcount(end) = mod(nt, ntperread);
end

nread = length(tstart);

fprintf('  Tranfering data:           \n');

for ii = 1:nread
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d of %3d\n', ii, nread);

    alldata = fread(fid, npertime.*tcount(ii), 'double');
    alldata = reshape(alldata, npertime, tcount(ii));
    alldata = mat2cell(alldata, sz, tcount(ii));

    for iv = 1:nvar

        if sz(iv) == 1
            start = tstart(ii);
            count = tcount(ii);
        else
            start = [0 tstart(ii)];
            count = [sz(iv) tcount(ii)];
        end
        
        netcdf.putVar(ncid, vid(iv), start, count, alldata{iv});

    end
    
    
end


