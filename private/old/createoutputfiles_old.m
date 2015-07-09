function ncvid = createoutputfiles(folder, z, zp, sd, ed, md, vartable, nens, newflag)

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

Dim(4).Name = 'ensemble';
Dim(4).Length = nens;
Dim(4).Unlimited = false;

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

if newflag

    if exist(folder, 'dir')
        warning('Output folder already exists; existing files will not be explicitly overwritten and may cause conflicts');
    else
        mkdir(folder);
    end

    dimfile = fullfile(folder, 'dimensions.nc');

    ncwriteschema(dimfile, C);

    ncwrite(dimfile, 'depth', z);
    ncwrite(dimfile, 'depth_edge', zp);
    ncwrite(dimfile, 'startdate', sd);
    ncwrite(dimfile, 'enddate', ed);
    ncwrite(dimfile, 'middate', md);
end

% Data variables

nvar = size(vartable, 1);

ncvid = zeros(nvar,2);

for iv = 1:nvar
    
    A = struct('Name', '/', 'Format', 'netcdf4');
    
    A.Variables.Name = vartable{iv,2};
    
    if isequal(size(vartable{iv,1}), [nz 1])
        A.Variables.Dimensions = Dim([1 3 4]);
    elseif isequal(size(vartable{iv,1}), [nzp 1])
        A.Variables.Dimensions = Dim([2 3 4]);
    elseif isequal(size(vartable{iv,1}), [1 1])
        A.Variables.Dimensions = Dim([3 4]);
    end
    
    A.Variables.Attributes = struct('Name', {'long_name', 'unit'}, 'Value', vartable(iv,3:4));
    
    A.Variables.Datatype = 'double';
    
    file = fullfile(folder, [vartable{iv,2} '.nc']);
    
    if newflag
        ncwriteschema(file, A);
    end
    
    ncvid(iv,1) = netcdf.open(file, 'WRITE');
    ncvid(iv,2) = netcdf.inqVarID(ncvid(iv,1), vartable{iv,2});
    
end

