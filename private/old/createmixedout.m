function createmixedout(filename, z, zp, sd, ed, md, vartable)

%-------------------
% Create file schema
%-------------------

A.Name = '/';
A.Format = 'netcdf4';

% Add dimensions

nz  = length(z);
nzp = length(zp);
nt  = size(sd,1);

A.Dimensions(1).Name = 'depth';
A.Dimensions(1).Length = nz;
A.Dimensions(1).Unlimited = false;

A.Dimensions(2).Name = 'depth_edge';
A.Dimensions(2).Length = nzp;
A.Dimensions(2).Unlimited = false;

A.Dimensions(3).Name = 'time';
A.Dimensions(3).Length = nt;
A.Dimensions(3).Unlimited = true;

% Add coordinate variables
    
A.Variables(1).Name = 'depth';
A.Variables(1).dshort = {'depth'};
A.Variables(1).ashort = {'unit', 'm', 'positive', 'up'};
A.Variables(1).Datatype = 'double';

A.Variables(2).Name = 'depth_edge';
A.Variables(2).dshort = {'depth_edge'};
A.Variables(2).ashort = {'unit', 'm', 'positive', 'up'};
A.Variables(2).Datatype = 'double';

A.Variables(3).Name = 'startdate';
A.Variables(3).dshort = {'time'};
A.Variables(3).ashort = {'unit', 'days since Jan-1-0000'};
A.Variables(3).Datatype = 'double';

A.Variables(4).Name = 'middate';
A.Variables(4).dshort = {'time'};
A.Variables(4).ashort = {'unit', 'days since Jan-1-0000'};
A.Variables(4).Datatype = 'double';

A.Variables(5).Name = 'enddate';
A.Variables(5).dshort = {'time'};
A.Variables(5).ashort = {'unit', 'days since Jan-1-0000'};
A.Variables(5).Datatype = 'double';

% Add data variables

for iv = 1:size(vartable,1)
    
    A.Variables(iv+5).Name = vartable{iv,2};
    
    if isequal(size(vartable{iv,1}), [nz 1])
        A.Variables(iv+5).dshort = {'depth', 'time'};
    elseif isequal(size(vartable{iv,1}), [nzp 1])
        A.Variables(iv+5).dshort = {'depth_edge', 'time'};
    elseif isequal(size(vartable{iv,1}), [1 1])
        A.Variables(iv+5).dshort = {'time'};
    end
    
    A.Variables(iv+5).ashort = {'long_name', vartable{iv,3}, 'unit', vartable{iv,4}};
    
    A.Variables(iv+5).Datatype = 'double';
end

% Add full info for attributes and dimensions

for iv = 1:length(A.Variables)
    [tf, loc] = ismember(A.Variables(iv).dshort, {A.Dimensions.Name});
    A.Variables(iv).Dimensions = A.Dimensions(loc);   
    
    A.Variables(iv).Attributes = cell2struct(reshape(A.Variables(iv).ashort, 2, []), {'Name', 'Value'});
end

A.Variables = rmfield(A.Variables, {'dshort', 'ashort'});

%-------------------
% Write to file
%-------------------

ncwriteschema(filename, A);

%-------------------
% Add coordinate 
% variable data
%-------------------

ncwrite(filename, 'depth', z);
ncwrite(filename, 'depth_edge', zp);
ncwrite(filename, 'startdate', sd);
ncwrite(filename, 'enddate', ed);
ncwrite(filename, 'middate', md);

% Test: fill data with nans?

% B = ncinfo(filename);
% 
% for iv = 6:length(B.Variables)
%     if isscalar(B.Variables(iv).Size)
%         sz = [B.Variables(iv).Size 1];
%     else
%         sz = B.Variables(iv).Size;
%     end
%     ncwrite(filename, B.Variables(iv).Name, nan(sz));
% end









