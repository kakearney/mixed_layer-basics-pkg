function [A, str] = initinterpdata(type, data, Grd)
%INITINTERPDATA Initialize datasets for later interpolation
%
% [A, str] = initinterpdata(type, data, Grd)
%
% This function sets up forcing datasets for use for later interpolation.
% It verifies that data is in the proper format (vs depth, time, or both),
% repeats climatological data if necessary, and extends data to the edges
% of the model temporal and spatial grids if necessary (using
% nearest-neighbor extrapolation).
%
% Input variables:
%
%   type:   One of the following strings:
%           'time':     data describes one or more variables vs time
%           'depth':    data describes one or more variables vs depth
%           'both':     data describes a single variable vs both time and
%                       depth
%
%   data:   Input data
%           'time':     nt x (6+nvar) array.  Columns 1-6 hold date vectors
%                       corresponding to each row of data, remaining
%                       column(s) hold variable data.  If the simulation
%                       time spans multiple years and the data only spans a
%                       single year, the data will be treated as a
%                       climatology and repeated for all simulation years
%           'depth':    nz x (1+nvar) array.  Column 1 holds depth values
%                       corresponding to each row of data (m, negative
%                       down), remaining column(s) hold variable data
%           'both':     (nt+1) x (nz+6) array.  Columns 1-6 hold date
%                       vectors corresponding to each row of data, row 1
%                       holds depth values corresponding to each column of
%                       data (cells (1,1:6) are ignored), remaining cells
%                       hold variable data.  If the simulation time spans
%                       multiple years and the data only spans a single
%                       year, the data will be treated as a climatology and
%                       repeated for all simulation years.
%
%   Grd:    Struct holding spatial and temporal grid data for mixed_layer
%           simulation
%
% Output variables:
%
%   A:      1 x 1 structure with the following fields:
%
%           t:      nt x 1 vector, time (seconds from simulation start
%                   time) ('time' or 'both' only) 
%
%           z:      nz x 1, depth (m, negative down) ('depth or 'both'
%                   only)
%
%           data:   nt x nvar ('time'), nz x nvar ('depth'), or nz x nt
%                   ('both') array of variable values 
%
%   str:    string, describing time and/or depth interval where
%           extrapolation was necessary

% Copyright 2009 Kelly Kearney

%--------------------------
% Setup
%--------------------------

t0 = datenum(Grd.start_date); 

tgap = zeros(0,2);
zgap = zeros(0,2);


switch type
    
    %--------------------------
    % Format time data
    %--------------------------
    
    case 'time'     
        
        % Separate data
        
        A.data = data(:,7:end);
        
        % Check for climatology
        
        dv = data(:, 1:6);
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dv, A.data, Grd.start_date(1), Grd.end_date(1), 2);
        end
        
        A.t = (datenum(dv) - t0) * 86400;
        
        % Extrapolate to edges of time grid
        
        if A.t(1) > 0
            tgap = [tgap; 0 A.t(1)];
            A.t = [0; A.t];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.t(end) < Grd.tmax
            tgap = [tgap; A.t(end) Grd.tmax];
            A.t = [A.t; Grd.tmax];
            A.data = [A.data; A.data(end,:)];
        end
        
    %--------------------------
    % Format depth data
    %--------------------------
        
    case 'depth'
        
        % Separate data
        
        A.z = data(:,1);
        A.data = data(:,2:end);
        
        % Extrapolate to edges of spatial grid
        
        if A.z(1) < 0
            zgap = [zgap; 0 A.z(1)];
            A.z = [0; A.z];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.z(end) > Grd.zp(end)
            zgap = [zgap; A.z(end) Grd.zp(end)];
            A.z = [A.z; Grd.zp(end)];
            A.data = [A.data; A.data(end,:)];
        end
        
    %--------------------------
    % Format time x depth data
    %--------------------------
    
    case 'both'
            
        % Separate data
            
        A.z = data(1,7:end)';
        A.data = data(2:end,7:end)';
        
        % Check for climatology
        
        dv = data(2:end, 1:6);
        if Grd.start_date(1) < Grd.end_date(1) && length(unique(dv(:,1))) == 1
            [dv, A.data] = repeatclimatology(dv, A.data, Grd.start_date(1), Grd.end_date(1), 1);
        end
        
        A.t = (datenum(dv) - t0) * 86400;
        
        % Extrapolate to edges of time and space grid
        
        if A.t(1) > 0
            tgap = [tgap; 0 A.t(1)];
            A.t = [0; A.t];
            A.data = [A.data(:,1) A.data];
        end
        
        if A.t(end) < Grd.tmax
            tgap = [tgap; A.t(end) Grd.tmax];
            A.t = [A.t; Grd.tmax];
            A.data = [A.data  A.data(:,end)];
        end
        
        if A.z(1) < 0
            zgap = [zgap; 0 A.z(1)];
            A.z = [0; A.z];
            A.data = [A.data(1,:); A.data];
        end
        
        if A.z(end) > Grd.zp(end)
            zgap = [zgap; A.z(end) Grd.zp(end)];
            A.z = [A.z; Grd.zp(end)];
            A.data = [A.data; A.data(end,:)];
        end

end

%--------------------------
% Extrapolation description
%--------------------------

str = extrapstr(tgap, t0, zgap);


%--------------------------
% Repeat climatological 
% data over simulation 
% period
%--------------------------

function [dvc, datac] = repeatclimatology(dv, data, syear, eyear, flag)

% if size(data,2) == size(dv,1) && size(data,1) ~= size(dv,1)
%     flag = 1;
% elseif size(data,1) == size(dv,1)
%     flag = 2;
% else
%     error('Can''t figure out data format');
% end

yrs = syear:eyear;
nyr = length(yrs);

if flag == 1    % depth x time data
    nperyear = size(data,2);
    dvc = repmat(dv, nyr, 1);
    dvc(:,1) = kron(yrs', ones(nperyear,1));
    datac = repmat(data, 1, nyr);
elseif flag == 2 % time x nvar data
    nperyear = size(data,1);
    dvc = repmat(dv, nyr, 1);
    dvc(:,1) = kron(yrs', ones(nperyear,1));
    datac = repmat(data, nyr, 1);
end

%--------------------------
% Create message telling 
% where extrapolation was
% needed
%--------------------------

function str = extrapstr(tgap, t0, zgap)
fmt = 'mm/dd/yyyy HH:MM';
dn = (tgap./86400) + t0;
ds = arrayfun(@(x) datestr(x,fmt), dn, 'uni', 0);

zstr = arrayfun(@(x) num2str(x, '%g m'), zgap, 'uni', 0);

str = [ds; zstr]';

str = sprintf('%s-%s, ', str{:});
str = str(1:end-2);







