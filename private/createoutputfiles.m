function [ncid, vid] = createoutputfiles(folder, z, zp, sd, ed, md, vartable, iens)

% Add dimensions

nz  = length(z);
nzp = length(zp);
nt  = size(sd,1);

Dim(1).Name = 'depth';
Dim(1).Length = nz;
Dim(1).Unlimited = false;

Dim(2).Name = 'depth_edge';
Dim(2).Length = nzp;
Dim(2).Unlimited = false;

Dim(3).Name = 'time';
Dim(3).Length = nt;
Dim(3).Unlimited = true;

% Dim(4).Name = 'ensemble';
% Dim(4).Length = nens;
% Dim(4).Unlimited = false;

% Coordinate variable file

C.Name = '/';
C.Format = 'netcdf4';

C.Variables(1).Name = 'depth';
C.Variables(1).Dimensions = Dim(1);
C.Variables(1).Attributes = struct('Name', {'unit', 'positive'}, 'Value', {'m', 'up'});
C.Variables(1).Datatype = 'double';

C.Variables(2).Name = 'depth_edge';
C.Variables(2).Dimensions = Dim(2);
C.Variables(2).Attributes = struct('Name', {'unit', 'positive'}, 'Value', {'m', 'up'});
C.Variables(2).Datatype = 'double';

C.Variables(3).Name = 'startdate';
C.Variables(3).Dimensions = Dim(3);
C.Variables(3).Attributes = struct('Name', {'unit'}, 'Value', {'days since Jan-1-0000'});
C.Variables(3).Datatype = 'double';

C.Variables(4).Name = 'middate';
C.Variables(4).Dimensions = Dim(3);
C.Variables(4).Attributes = struct('Name', {'unit'}, 'Value', {'days since Jan-1-0000'});
C.Variables(4).Datatype = 'double';

C.Variables(5).Name = 'enddate';
C.Variables(5).Dimensions = Dim(3);
C.Variables(5).Attributes = struct('Name', {'unit'}, 'Value', {'days since Jan-1-0000'});
C.Variables(5).Datatype = 'double';

% Create folder

if exist(folder, 'dir')
%     warning('Output folder already exists; existing files will not be explicitly overwritten and may cause conflicts');
else
    mkdir(folder);
end

dimfile = fullfile(folder, 'dimensions.nc');
if ~exist(dimfile, 'file')

    try
        ncwriteschema(dimfile, C);
    catch
        pause(5); % In case of weird parallel bug, pause and try again
        W = warning('off', 'all'); % b/c nc.createvariable is going to throw a bazillion warnings about variables existing
        ncwriteschema(dimfile, C);
        warning(W);
    end

    ncwrite(dimfile, 'depth', z);
    ncwrite(dimfile, 'depth_edge', zp);
    ncwrite(dimfile, 'startdate', sd);
    ncwrite(dimfile, 'enddate', ed);
    ncwrite(dimfile, 'middate', md);
end

% Data variables

nvar = size(vartable, 1);

ncvid = zeros(nvar,2);

A = struct('Name', '/', 'Format', 'netcdf4');

for iv = 1:nvar
    
    A.Variables(iv).Name = vartable{iv,2};
    
    if isequal(size(vartable{iv,1}), [nz 1])
        A.Variables(iv).Dimensions = Dim([1 3]);
    elseif isequal(size(vartable{iv,1}), [nzp 1])
        A.Variables(iv).Dimensions = Dim([2 3]);
    elseif isequal(size(vartable{iv,1}), [1 1])
        A.Variables(iv).Dimensions = Dim([3]);
    end
    
    A.Variables(iv).Attributes = struct('Name', {'long_name', 'unit'}, 'Value', vartable(iv,3:4));
    
    A.Variables(iv).Datatype = 'double';
    A.Variables(iv).FillValue = -9999;
    
end

file = fullfile(folder, sprintf('sim%04d.nc', iens)); %[vartable{iv,2} '.nc']);    
if exist(file, 'file')
    warning('WCE:fileexists', 'Specified outputfile/iens combo already exists; attempting to overwrite file');
    W = warning('off', 'all'); % nc.createvariable thing again
    ncwriteschema(file, A);
    warning(W);
else
    ncwriteschema(file, A);
end  
    
ncid = netcdf.open(file, 'WRITE');
vid = zeros(nvar,1);
for iv = 1:nvar
    vid(iv) = netcdf.inqVarID(ncid, vartable{iv,2});
end

