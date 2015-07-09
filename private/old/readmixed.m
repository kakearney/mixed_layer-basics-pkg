function Data = readmixed(file, varargin)
%READMIXED Read mixed_layer output into a structure
%
% Data = readmixed(file)
% Data = readmixed(file, 'var1', 'var2', ...)
% Data = readmixed(file, Opt, 'var1', 'var2', ...)
% readmixed(file, 'list')
%
% Input variables:
% 
%   file:       name of output netcdf file (with or without .nc extension)
%
%   var#:       name of variables to be read.  If not included, all
%               variables will be read.  If the string 'dimensions' is
%               included, this corresponds to all the dimensional
%               variables: 'depth', 'depth_edge', 'startdate', 'middate',
%               'enddate'.
%
%   Opt:        1 x 1 one structure holding options for reading and
%               averaging data:
%
%               tstart:     date number corresponding to earliest value to
%                           read in time dimension.  If NaN, will read from
%                           beginning. [NaN]
%
%               tend:       date number corresponding to latest value to
%                           read in time dimension.  If NaN, data will be
%                           read to end of time dimension. [NaN]
%
%               tstride:    sampling interval (stride) along time
%                           dimension. [1]
%
%               tbin:       Bin edges, as in histc.  Data will be averaged
%                           over each bin.  Sampling based on above
%                           parameters will be done prior to averaging.
%
%               monthavg:   True indicates to average data over each month
%                           (a shortcut to listing months as tbin).
%                           Sampling based on above paramters will be done
%                           prior to averaging. [false]
%
%               tmid:       Values of the 'middate' variable in the file.
%                           If you're doing multiple reads of a large file
%                           using tstart or tend, it can speed things up if
%                           readmixed doesn't need to read in this variable
%                           each time to analyze where to start and stop.
%
%               surfonly:   Only read the surface depth values, i.e. first
%                           cell in the depth dimension [false] 
%
%   'list':     When called with just the file name and this string, this
%               will simply list the names and dimensions of all variables
%               in the file.
%
% Output variables:
%
%   Data:       1 x 1 structure holding data from specified file.  Field
%               names correspond to variable names. 

% Copyright 2009 Kelly Kearney

%-----------------------------
% Check input
%-----------------------------

[pth,fl,ext] = fileparts(file);
if isempty(ext)
    file = [file '.nc'];
end

% if ~strcmp(file(end-2:end), '.nc')
%     file = [file '.nc'];
% end
if ~exist(file, 'file')
    error('Results file not found');
end    

%-----------------------------
% List variables
%-----------------------------

Info = nc_info(file);

% usenew = ~verLessThan('matlab', '7.10.0');
usenew = false;

% Some versions of nc_info use DataSet, some Dataset

fld = fieldnames(Info);
isds = strcmpi('dataset', fld);
dsfld = fld{isds};

varnames = {Info.(dsfld).Name};

if nargin == 2 && strcmp(varargin{1}, 'list')
    
    dims = {Info.(dsfld).Dimension};
    
    isz = strcmp(varnames, 'depth');
    nz = Info.(dsfld)(isz).Size;
    iszp = strcmp(varnames, 'depth_edge');
    nzp = Info.(dsfld)(iszp).Size;
    ist = strcmp(varnames, 'startdate');
    nt = Info.(dsfld)(ist).Size;
    
    dimstr = cellfun(@(x) sprintf('%s x ', x{:}), dims, 'uni', 0);
    dimstr = cellfun(@(x) x(1:end-3), dimstr, 'uni', 0);
    
    dimstr{isz} = sprintf('%d x 1', nz);
    dimstr{iszp} = sprintf('%d x 1', nzp);
    dimstr{ist} = sprintf('%d x 1', nt);
   
    data = [{'Variable', 'Dimensions'; '--------', '----------'}; [varnames' dimstr']]';
    
    if nargout == 1
        Data = data(:,3:end)';
    else
        fprintf('%-10s  %-20s\n', data{:});
    end
else
    
    % Parse options and variable names
    
    Opt.tstride = 1;
    Opt.tstart = NaN;
    Opt.tend = NaN;
    Opt.monthavg = false;
    Opt.yearavg = false;
    Opt.tbin = [];
    Opt.tmid = [];
    Opt.surfonly = false;
    
    if nargin == 1
        vars = varnames;
    else
        isopt = cellfun(@isstruct, varargin);
        if any(isopt)
            if sum(isopt) > 1
                error('Only one options structure allowed');
            end
            pv(1,:) = fieldnames(varargin{isopt});
            pv(2,:) = struct2cell(varargin{isopt});
            Opt = parsepv(Opt, pv(:));
%             S = warning('off');
%             Opt = mergestruct(Opt, varargin{isopt});
%             warning(S);
            vars = varargin(~isopt);
        else
            vars = varargin;
        end
    end

%         isscalar(varargin{1}) && isnumeric(varargin{1})
%         timestride = varargin{1};
%         if nargin > 2
%             vars = varargin(2:end);
%         else
%             vars = varnames;
%         end
%     else
%         timestride = 1;
%         vars = varargin;
%     end
    
    % Expand dimension variable names if necessary

    isdimshortcut = strcmp(vars, 'dimensions');
    if any(isdimshortcut)
        vars = vars(~isdimshortcut);
        dimnames = {'depth', 'depth_edge', 'startdate', 'middate', 'enddate'};
        [tf,loc] = ismember(vars, dimnames);
        vars = vars(~tf);
        vars = [dimnames vars];
    end
     
    % Check that all variables exist in file 
    
    [tf, loc] = ismember(vars, varnames);
    if ~all(tf)
        varstr = sprintf('%s, ', vars{~tf});
        error('Did not find variable(s) in file: %s', varstr(1:end-2));
    end
    
    % Figure out start, count, and stride values for time dimension
    
    if usenew
        ncid = netcdf.open(file, 'NC_NOWRITE');
    end
    
    needtime = ~isnan(Opt.tstart) || ~isnan(Opt.tend) || Opt.monthavg || Opt.yearavg;
    if needtime
        if isempty(Opt.tmid)
            if usenew
                tmidid = netcdf.inqVar(ncid, 'middate');
                tmid = netcdf.getVar(ncid, tmidid);
            else
                tmid = nc_varget(file, 'middate');
            end
            nt = length(tmid);
        else
            tmid = Opt.tmid;
            istime = strcmp({Info.Dimension.Name}, 'time');
            nt = Info.Dimension(istime).Length;
            if length(tmid) ~= nt
                error('tmid length doesn''t match the middate variable in file');
            end       
        end
            
    else
        istime = strcmp({Info.Dimension.Name}, 'time');
        nt = Info.Dimension(istime).Length;
    end

    if isnan(Opt.tstart)
        tstart = 0;
    else
        tstart= find(tmid > Opt.tstart, 1, 'first') - 1;
    end
        
    if Opt.tstride == 1 && isnan(Opt.tend)
        tcount = -1;
    else
        if isnan(Opt.tend)
            tend = nt - 1;
        else
            tend = find(tmid < Opt.tend, 1, 'last') - 1;
        end
        npt = tend - tstart + 1;
        tcount = floor(npt/Opt.tstride);
    end
            
    % Read variables from file
    
    if Opt.yearavg || Opt.monthavg || ~isempty(Opt.tbin)
        if usenew
            tmidid = netcdf.inqVar(ncid, 'middate');
            tmid = netcdf.getVar(ncid, tmidid, tstart, tcount, Opt.tstride);
        else
            tmid = nc_varget(file, 'middate', tstart, tcount, Opt.tstride);
        end
        dv = datevec(tmid);
    end
    
    for iv = 1:length(vars)
        
        vardims = Info.(dsfld)(loc(iv)).Dimension;
        istime = strcmp(vardims, 'time');
        isdepth = strcmp(vardims, 'depth');
        nt = Info.(dsfld)(loc(iv)).Size(istime);
        
        start = zeros(size(vardims));
        count = ones(size(vardims)) * -1;
        stride = ones(size(vardims));
        
        start(istime) = tstart;
        count(istime) = tcount;
        stride(istime) = Opt.tstride;
        
        if Opt.surfonly
            start(isdepth) = 0;
            count(isdepth) = 1;
        end
           
        if usenew
            vid = netcdf.inqVarID(ncid,vars{iv});
            Data.(vars{iv}) = netcdf.getVar(ncid, vid, start, count, stride);
        else
            Data.(vars{iv}) = nc_varget(file, vars{iv}, start, count, stride);
        end
       
        
        % Monthly averages, if requested
        
        if Opt.monthavg
            if any(istime)
                if strcmp(vars{iv}, 'startdate')
                    [blah, Data.(vars{iv})] = consolidator(dv(:,1:2), Data.(vars{iv}), @min);
                elseif strcmp(vars{iv}, 'enddate')
                    [blah, Data.(vars{iv})] = consolidator(dv(:,1:2), Data.(vars{iv}), @max);
                else
                    [blah, Data.(vars{iv})] = consolidator(dv(:,1:2), Data.(vars{iv}), @nanmean);
                end
            end
        elseif ~isempty(Opt.tbin)
            if any(istime)
                [n,bin] = histc(tmid, Opt.tbin);
                if strcmp(vars{iv}, 'startdate')
                    [blah, Data.(vars{iv})] = consolidator(bin, Data.(vars{iv}), @min);
                elseif strcmp(vars{iv}, 'enddate')
                    [blah, Data.(vars{iv})] = consolidator(bin, Data.(vars{iv}), @max);
                else
                    [blah, Data.(vars{iv})] = consolidator(bin, Data.(vars{iv}), @nanmean);
                end
            end
        end
                
    end
    
    % Remove nulls

    fld = fieldnames(Data);
    for ifld = 1:length(fld)
        isnull = Data.(fld{ifld}) > 1e36;
        Data.(fld{ifld})(isnull) = NaN;
    end
    
end



    
    

