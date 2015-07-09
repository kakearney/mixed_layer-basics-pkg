function binary2netcdf(binfile, ncfile, sdate, edate, mdate, z, zp, vartable, ntperfile, verbose)
%BINARY2NETCDF Convert temporary binary output file to netcdf
%
% binary2netcdf(binfile, ncfile, sdate, edate, mdate, z, zp, vartable, ...
%               ntperread)
%
% This function creates the netcdf output file for mixed_layer.  It
% requires Matlab R2009a or later, which has netcdf support.
% binary2netcdf_snc can replace this function for older Matlab versions,
% but that version is currently absurdly slow. 
%
% Input variables:
%
%	binfile:    name of binary output file
%
%	ncfile:     name of netcdf file
%
%	sdate:      vector of datenumbers, time grid (beginning) 
%
%	edate:      vector of datenumbers, time grid (end)
%
%	mdate:      vector of datenumbers, time grid (middle)
%
%	z:          depth grid
%
%	zp:         depth grid, edges
%
%	vartable:   n x 3 cell array, all variables from simulation
%
%	ntperfile:  number of time steps to be written to each small
%               temporary file
%
%	verbose:    flag indicating whether to print progress to screen

% Copyright 2010 Kelly Kearney

sz = cellfun(@(x) size(x,1), vartable(:,1));
npertime = sum(sz);
nvar = length(sz);

% How many time steps were written to this file?

fid = fopen(binfile);
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

% Figure out file and read-write breakdown

timeidx = wreshape(1:nt, ntperfile, [], 'pad', NaN);
ntinfile = sum(~isnan(timeidx));
nfile = size(timeidx, 2);

ntperread = 30;
readidx = zeros(0,3); % col 1 = index of file, col 2 = number of times to read
for ifl = 1:nfile    
    nfull = floor(ntinfile(ifl)./ntperread);
    partial = mod(ntinfile(ifl), ntperread);
    temp = ones(nfull,1)*[ifl ntperread];
    if partial > 0
        temp = [temp; ifl partial];
    end
    readidx = [readidx; [temp (0:ntperread:ntinfile(ifl)-1)']];
end

% Create temporary files

if verbose
	fprintf('  Creating file:           \n');
end

for ifl = 1:nfile
    
	if verbose
    	fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d of %3d\n', ifl, nfile);
	end
    
    filename = regexprep(ncfile, '\.nc', sprintf('_%03d.tmp.nc', ifl));
    idx = timeidx(1:ntinfile(ifl),ifl);        
    [fileid(ifl), varid] = createfile(filename, z, zp, sdate(idx,:), edate(idx,:), mdate(idx,:), vartable);
end

% Read in data from binary file and write to netcdf

nread = size(readidx,1);

if verbose
	fprintf('  Writing data:           \n');
end

for ii = 1:nread
    
	if verbose
    	fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d of %3d\n', ii, nread);
	end

    alldata = fread(fid, npertime.*readidx(ii,2), 'double');
    alldata = reshape(alldata, npertime, readidx(ii,2));
    alldata = mat2cell(alldata, sz, readidx(ii,2));

    for iv = 1:nvar

        if sz(iv) == 1
            start = readidx(ii,3);
            count = readidx(ii,2);
        else
            start = [0      readidx(ii,3)];
            count = [sz(iv) readidx(ii,2)];
        end
        
        netcdf.putVar(fileid(readidx(ii,1)), varid(iv), start, count, alldata{iv});
    end
    
    
end

% Close access to temporary files and binary file

for ifl = 1:nfile
    netcdf.close(fileid(ifl));
end 

fclose(fid);

% Combine temporary files into one

if verbose
	fprintf('  Concatenating files\n');
end

if nfile > 1
    cmd = sprintf('ncrcat -O %s_[0-9][0-9][0-9].tmp.nc %s', ncfile(1:end-3), ncfile);
    [s,r] = system(cmd);
    if s
        error('Error concatenating files:%s', r);
    end
    for ifl = 1:nfile
        delete(sprintf('%s_%03d.tmp.nc', ncfile(1:end-3), ifl));
    end
else
    movefile(sprintf('%s_001.tmp.nc', ncfile(1:end-3)), ncfile);
end
    



function [ncid, vid] = createfile(filename, z, zp, sd, ed, md, vartable)
% Start file

mode = bitor(netcdf.getConstant('NC_CLOBBER'), netcdf.getConstant('NC_64BIT_OFFSET')); 
ncid = netcdf.create(filename, mode);

% Add dimensions

nz  = length(z);
nzp = length(zp);
nt  = size(sd,1);

zid  = netcdf.defDim(ncid, 'depth', nz);
zpid = netcdf.defDim(ncid, 'depth_edge', nzp);
tid  = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));

% Coordinate variables

cvid(1) = netcdf.defVar(ncid, 'depth',      'double', zid);
cvid(2) = netcdf.defVar(ncid, 'depth_edge', 'double', zpid);
cvid(3) = netcdf.defVar(ncid, 'startdate',  'double', tid);
cvid(4) = netcdf.defVar(ncid, 'enddate',    'double', tid);
cvid(5) = netcdf.defVar(ncid, 'middate',    'double', tid);

netcdf.putAtt(ncid, cvid(1), 'unit', 'm');
netcdf.putAtt(ncid, cvid(1), 'positive', 'up');
netcdf.putAtt(ncid, cvid(2), 'unit', 'm');
netcdf.putAtt(ncid, cvid(2), 'positive', 'up');
netcdf.putAtt(ncid, cvid(3), 'unit', 'days since Jan-1-0000');
netcdf.putAtt(ncid, cvid(4), 'unit', 'days since Jan-1-0000');
netcdf.putAtt(ncid, cvid(5), 'unit', 'days since Jan-1-0000');

netcdf.endDef(ncid);

netcdf.putVar(ncid, cvid(1), z);
netcdf.putVar(ncid, cvid(2), zp);
netcdf.putVar(ncid, cvid(3), 0, nt, sd);
netcdf.putVar(ncid, cvid(4), ed);
netcdf.putVar(ncid, cvid(5), md);

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






