 function varargout = archivemldata(Grd, Arch, outfile, it, data, verbose)
%ARCHIVEMLDATA Write mixed_layer results to file
%
% avg = archivemldata(Grd, Arch, outfile, it, data)
% [avg, fid] = archivemldata(Grd, Arch, outfile, it, data)
%
% New version should be compatible with all versions of Matlab, as long as
% snctools is installed and mexnc points to the proper backend (either
% newer internal stuff or the proper mex version).
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
%               or scalars.  Columns 2-4 are strings with the variables'
%               short names, long names, and units, respectively.
%
% Output variables:
%
%   avg:        nvar x 1 cell array, holding the newly-averaged values each
%               variable.  Averages are performed over each archiving
%               period.

% Copyright 2010 Kelly Kearney


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

% Write to temporary binary file

if it == 1
    tempfile = tempname;
    if ~isempty(Arch.tempdir)
        [blah, tempfile] = fileparts(tempfile);
        tempfile = fullfile(Arch.tempdir, tempfile);
    end
    mlmarker = 'mixed_layer temp file';
    save(tempfile, 'data', 'outfile', 'mlmarker'); % .mat file
    fid = fopen([tempfile '.bin'], 'w');
    if Arch.islast(1)
        fwrite(fid, cat(1, avg{:}), 'double');
        avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
    end
elseif Arch.islast(it)
    fwrite(Arch.fid, cat(1, avg{:}), 'double');
    
    % Zero out avg arrays
    
    avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
end

% Assign output

varargout{1} = avg;
if nargout > 1
    varargout{2} = fid;
    varargout{3} = tempfile;
end

% Convert from binary file to netcdf file

if it == Grd.nt
   
    fclose(Arch.fid);

	if verbose
    	fprintf('Saving simulation results...\n');
	end	
    
    sdate = Arch.dateedge(1:end-1);
    edate = Arch.dateedge(2:end);
    
    binary2netcdf([Arch.tempfile '.bin'], outfile, sdate, edate, ...
                  Arch.middate, Grd.z, Grd.zp, data, Arch.nper, verbose);
              
    if Arch.cleanup
       delete([Arch.tempfile '.bin']);
       delete([Arch.tempfile '.mat']);
    end
    
end




% Create new output file, if beginning of new file timestep

% if Arch.isnewfile(it)
%     
%     % Create new file
%     
%     timeidx = Arch.filedates(Arch.fileidx(it,1),:);
%     timeidx = timeidx(1):timeidx(2);
%     sdate = Arch.startdate(timeidx,:); %TODO add initial time step?
%     edate = Arch.enddate(timeidx,:);
%     mdate = Arch.middate(timeidx,:);
%     
%     createfile(outfile, Arch.fileidx(it), sdate, edate, mdate, data, Grd.z, Grd.zp, timeidx);
%     
% end

% Save to file, if end of averaging time step
    
% if Arch.islast(it)
%     
%     filename = regexprep(outfile, '\.nc', sprintf('_%03d.tmp.nc', Arch.fileidx(it)));
%     bin = Arch.fileidx(it,2);
%    
%         
%     % Write to file (snctools method, much faster on my computer)
%     
%     for ivar = 1:size(data,1)
%         if isscalar(avg{ivar})
%             start = bin-1;
%             count = 1;
%         elseif isvector(avg{ivar});
%             start = [bin-1 0];
%             count = [1 length(avg{ivar})];
%         end
%         nc_varput(filename, data{ivar,2}, avg{ivar}, start, count);
%     end
%     
%     % Zero out avg arrays
%     
%     avg = cellfun(@(x) zeros(size(x)), avg, 'uni', 0);
%     
% end
% 
% % Combine files if reached last time step
% 
% if it == Grd.nt
%     
%     nfile = size(Arch.filedates,1);
%     
%     if nfile > 1
%         cmd = sprintf('ncrcat -O %s_[0-9][0-9][0-9].tmp.nc %s', outfile(1:end-3), outfile);
%         [s,r] = system(cmd);
%         if s
%             error('Error concatenating files:%s', r);
%         end
%         for ifl = 1:nfile
%             delete(sprintf('%s_%03d.tmp.nc', outfile(1:end-3), ifl));
%         end
%     else
%         movefile(sprintf('%s_001.tmp.nc', outfile(1:end-3)), outfile);
%     end
%     
% end
    


% function createfile(outfile, ifl, sdate, edate, mdate, vartable, z, zp, times)
% 
% % Start file
% 
% filename = regexprep(outfile, '\.nc', sprintf('_%03d.tmp.nc', ifl));
% 
% mode = bitor(nc_clobber_mode, nc_64bit_offset_mode);
% nc_create_empty(filename, mode);
% 
% % Add dimensions
% 
% nz = length(z);
% nzp = length(zp);
% nt = size(sdate,1);
% 
% nc_add_dimension(filename, 'depth', nz);
% nc_add_dimension(filename, 'depth_edge', nzp);
% nc_add_dimension(filename, 'time', 0);
% 
% % Coordinate variables
% 
% Var.Name = 'depth';
% Var.Dimension = {'depth'};
% nc_addvar(filename, Var);
% nc_varput(filename, Var.Name, z);
% nc_attput(filename, Var.Name, 'unit', 'm');
% 
% Var.Name = 'depth_edge';
% Var.Dimension = {'depth_edge'};
% nc_addvar(filename, Var);
% nc_varput(filename, Var.Name, zp);
% nc_attput(filename, Var.Name, 'unit', 'm');
% 
% Var.Name = 'startdate';
% Var.Dimension = {'time'};
% nc_addvar(filename, Var);
% nc_varput(filename, Var.Name, datenum(sdate));
% nc_attput(filename, Var.Name, 'unit', 'serial date number');
% 
% Var.Name = 'enddate';
% Var.Dimension = {'time'};
% nc_addvar(filename, Var);
% nc_varput(filename, Var.Name, datenum(edate));
% nc_attput(filename, Var.Name, 'unit', 'serial date number');
% 
% Var.Name = 'middate';
% Var.Dimension = {'time'};
% nc_addvar(filename, Var);
% nc_varput(filename, Var.Name, datenum(mdate));
% nc_attput(filename, Var.Name, 'unit', 'serial date number');
% 
% % Data variables
% 
% for ivar = 1:size(vartable,1)
%     
%     Var.Name = vartable{ivar,2};
%     
%     if isequal(size(vartable{ivar,1}), [nz 1])
%         Var.Dimension = {'time', 'depth'};
%         start = [0 0];
%         count = [1 nz];
%     elseif isequal(size(vartable{ivar,1}), [nzp 1])
%         Var.Dimension = {'time', 'depth_edge'};
%         start = [0 0];
%         count = [1 nzp];
%     elseif isequal(size(vartable{ivar,1}), [1 1])
%         Var.Dimension = {'time'};
%         start = 0;
%         count = 1;
%     else
%         error('Archive variable has unexpected dimensions');
%     end
%     
%     nc_addvar(filename, Var);
%     nc_varput(filename, Var.Name, vartable{ivar,1}, start, count);
%     
%     nc_attput(filename, Var.Name, 'long_name', vartable{ivar,3});
%     nc_attput(filename, Var.Name, 'unit', vartable{ivar,4});
%     
% end

