function [Grd, Ht, Ts, Mmntm, Wnd, Arch] = initialize(In)
%INITIALIZE Initialize parameters for mixed_layer model
%
% This function initializes the values of most constant parameters used
% throughout the mixed-layer model.  It also preallocates and sets the
% initial conditions for many variables that will be modified throughout
% the model simulation.
%
% Input variables:
%
%   In:     structure holding user-supplied input variables
%
% Output variables:
%
%   Grd:    structure holding temporal and spatial grid parameters:
%
%           z:          nz x 1 array, depth coordinate at the center of
%                       each grid cell (m) 
%
%           zp:         (nz+1) x 1 array, depth coordinate at the edges of
%                       each grid cell (m) 
%
%           nz:         number of vertical levels 
%
%           tmax:       simulation length (s)
%
%           nt:         total number of internal time iterations
%
%           time:       1 x nt array, time elapsed from model start time to
%                       the beginning of each time interval (s) 
%
%           start_date: 1 x 6 array, date vector for simulation start date.
%                       This will always be Jan 1 of the specified start
%                       year.
%
%           end_date:   1 x 6 array, date vector for simulation end date.
%                       This will always be Dec 31 of the specified end
%                       year.
%
%   Ht:     structure holding variables related to heat forcing
%
%           t:          nth x 1 array, time corresponding to Ht.data (s
%                       from simulation start)
%
%           data:       nth x 3 array, heat forcing data.  Column 1 =
%                       incoming solar radiation (W m^-2), Column 2 =  air
%                       temperature (deg C), Column 3 = dewpoint
%                       temperature (deg C)
%
%           Qo:         nt x 1 array, estimate of clear sky irrandiance,
%                       based on the "Smithsonian Formula" from Seckel and
%                       Beaudry, as reported in Reed, 1977, JPO, 7, pp.
%                       482-485.  It is good for latitudes between 20S and
%                       60N.
%
%           meanQi:     nt x 1 array, mean observed daily irradiance (W
%                       m^-2)
%
%   Ts:     structure holding variables related to temperature and salinity
%
%           T:          nz x 1 array, temperature profile (deg C)
%
%           S:          nz x 1 array, salinity profile (psu)
%
%           Sig:        nz x 1 array, density profile.  Note that water is
%                       currently treated as incompressible (kg m^-3)
%
%           Srelax:     1 x 1 structure of data for salt relaxation
%                       interpolation.  Not included if no salt relaxation
%                       data was provided.
%                       t:      nts x 1 array, times corresponding to
%                               columns of data (s from sim start time) 
%                       z:      nzs x 1 array, depths corresponding to the
%                               rows of data (m, neg down) 
%                       data:   nzs x nts array, salt relaxation profiles
%                               (psu) 
%
%           Trelax:     1 x 1 structure of data for temperature relaxation
%                       interpolation.  Not included if no temperature
%                       relaxation data was provided.
%                       t:      ntt x 1 array, times corresponding to
%                               columns of data (s from sim start time) 
%                       z:      nzt x 1 array, depths corresponding to the
%                               rows of data (m, neg down) 
%                       data:   nzt x ntt array, temperature relaxation
%                               profiles (deg C)  
%
%   Mmntm:  structure holding variables related to momentum, mixing, and
%           turbulence 
%
%           Kh:         (nz+1) x 1 array, tracer mixing coefficient (m^2
%                       s^-1) 
%
%           small:      arbitrarily small value, used as proxy for 0 to
%                       prevent things from blowing up 
%
%           U:          nz x 1 array, east-west current velocity (m/s)
%
%           V:          nz x 1 array, north-south current velocity (m/s)
%
%           q2:         (nz+1) x 1 array, twice the turbulent kinetic
%                       energy (m^2 s^-2)
%
%           q2l:        (nz+1) x 1 array, turbulent kinetic energy * length
%                       scale term (m^3 s^-2)  
%
%           len:        (nz+1) x 1 array, turbulence length scale (m)
%
%           gh:         (nz+1) x 1 array, Richardson number (no units)
%
%           sh:         (nz+1) x 1 array, MY 2.5 intermediate quantity, a
%                       function of the Richardson number
%
%           sm:         (nz+1) x 1 array, MY 2.5 intermediate quantity, a
%                       function of the Richardson number
%
%           Km:         (nz+1) x 1 array, vertical kinematic viscosity,
%                       i.e. momentum mixing coefficient (m^2 s^-1) 
%
%           Kq:         (nz+1) x 1 array, turbulence mixing coefficient
%                       (m^2 s^-1)
%
%           boygr:      (nz+1) x 1 array, buoyancy generation term
%
%           shear:      (nz+1) x 1 array, shear term
%
%           lc_q2:      (nz+1) x 1 array, dissipation constant for
%                       turbulent kinetic energy (no units)
%
%           lc_q2l:     (nz+1) x 1 array, dissipation constant for
%                       kinetic-energy-length-scale (no units)
%
%           kmol:		Background diffusivity (m^2 s^-1)	
%
%           cor:        coriolis forcing based on latitude (s^-1)
%
%           pgx:        acceleration due to a pressure gradient in the EW
%                       direction (m s^-2)
%
%           pgy:        accelaration due to a pressure gradient in the NS
%                       direction (m s^-2)
%
%           kappa:      parameter related to bottom friction, von Karman's
%                       constant (no units)
%
%           z0b:        roughness parameter, related to bottom friction (m)
% 
%           Cbot:       bottom friction coefficient (POM parameterization)
%                       (no units)
%
%   Wnd:    structure holding variables related to wind forcing
%
%           t:          ntw x 1 array, times corresponding to wind data (s
%                       from sim start time)
%
%           data:       ntw x 3 array of wind forcing data. Column 1 =
%                       surface wind stress in east-west (u) direction (N
%                       m^-2), column 2 = surface wind stress in the
%                       north-south (v) direction (N m^-2), column 3 = wind
%                       speed at 10 m above sea level (m/s)
%
%   Arch:   structure holding variables related to archive (i.e. output).
%           The archiving period refers to the In.tarch input.
%
%           startdate:  nbin x 6 array, date vectors corresponding to
%                       start of each archiving time step
%
%           enddate:    nbin x 6 array, date vectors corresponding to
%                       end of each archiving time step
%
%           middate:    nbin x 6 array, date vectors corresponding to
%                       middle of each archiving time step
%
%           fraction:   1 x nt array, fraction that each model time step
%                       contributes to its archiving time step
%
%           islast:     1 x nt logical array, true if time step
%                       corresponds to the end of an archiving time step
%
%           bin:        1 x nt array, index of archiving time step to which
%                       each model time step corresponds
%
%           nbin:       number of archiving time steps
%
%           fileidx:    nt x 2 array, column one holds the index of the
%                       temporary output file to which results will be
%                       written for each time step, column 2 tells the
%                       index of that time step within the file
%
%           isnewfile:  nt x 1 array, true if model time step is the first
%                       to be written to a new temporary output file
% 
%           endtime:    n x 1, end time of each archiving period (s)
%
%           filedates:  nfile x 2 array, indices of first and last archive
%                       step included in each temporary output file
%           
%
% Charlie Stock
% cstock@alum.mit.edu
%
% modified by Kelly Kearney

%--------------------------
% Calculate the grid and 
% time step properties
%--------------------------

if isscalar(In.dz)
    Grd.z = [-In.dz/2:-In.dz:In.zbot+In.dz/2]'; 
    Grd.zp = [0:-In.dz:In.zbot]';            
    Grd.nz = length(Grd.z);
else
    In.dz = In.dz(:);
    Grd.zp = -[0; cumsum(In.dz)];
    Grd.z = (Grd.zp(1:end-1) + Grd.zp(2:end))./2;
    Grd.nz = length(Grd.z);
end

if isequal(size(In.syear), [1 6])
    Grd.start_date = In.syear;
else
    Grd.start_date = [In.syear 1 1 0 0 0];      % Start at midnight morning Jan 1
end

if isequal(size(In.eyear), [1 6])
    Grd.end_date = In.eyear;
else
    Grd.end_date   = [In.eyear 12 31 24 0 0];   % End midnight night Dec. 31
end

dnstart = datenum(Grd.start_date);
dnend = datenum(Grd.end_date);

Grd.tmax = (dnend - dnstart)*86400;
Grd.nt = floor(Grd.tmax/In.dt);

Grd.time = (0:Grd.nt-1)*In.dt;

% For screen print counter: print progress at the end of each day

daydiff = diff(floor(Grd.time/86400 + datenum(Grd.start_date)));
Grd.newday = [true logical(daydiff)];

% Hotstart save index

if ~isempty(In.hotstartdn)
    Grd.savehot = find(Grd.time < (In.hotstartdn - dnstart)*86400, 1, 'last');
else
    Grd.savehot = NaN;
end

%--------------------------
% Interpolate the wind and 
% heat input onto the model 
% time steps
%--------------------------

[Ht, str] = initinterpdata('time', In.heat_input, Grd);   % Qi, airtemp, dew point
if ~isempty(str) && In.verbose
    fprintf('  Missing heat data: %s\n', str);
end

[WndSpd, str] = initinterpdata('time', In.wind_input, Grd); % uwspeed, vwspeed
if ~isempty(str) && In.verbose
    fprintf('  Missing wind data: %s\n', str);
end


%--------------------------
% Initialize temperature, 
% salinity and density 
% profiles and, if using 
% the Mellor-Yamada 
% turbulence closure, set 
% the velocity to 0. Also, 
% set Boundary conditions 
% for temperature and 
% salinity                                                       %
%--------------------------

TsIntp = initinterpdata('depth', In.ts_input, Grd);

ts = interp1(TsIntp.z, TsIntp.data, Grd.z);
Ts.T = ts(:,1);
Ts.S = ts(:,2);

if any(isnan(Ts.T(:)))
    error('NaN found in temperature data');
end
if any(isnan(Ts.S(:)))
    error('NaN found in salinity data');
end

Ts.Sig = sw_dens0(Ts.S, Ts.T);


%--------------------------
% Set up temperature and
% salinity relaxation
%--------------------------

% Salinity

if In.hassrelax
    
    [Ts.Srelax, str] = initinterpdata('both', In.srelax, Grd);
    
    if ~isempty(str) && In.verbose
        fprintf('  Missing salt relaxation data: %s\n', str);
    end
end

% Temperature

if In.hastrelax
    [Ts.Trelax, str] = initinterpdata('both', In.trelax, Grd);
    
    if ~isempty(str) && In.verbose
        fprintf('  Missing temperature relaxation data: %s\n', str);
    end
end


%--------------------------
% T/S/tracer upwelling or
% downwelling
%--------------------------

if In.hasw
    
    [Ts.Upwell, str] = initinterpdata('both', In.tracerw, Grd);
    
    if ~isempty(str) && In.verbose
        fprintf('  Missing up/downwell data: %s\n', str);
    end
   
end

%--------------------------
% Initialize mixing and 
% turbulence variables
%--------------------------

Mmntm.Kh = zeros(Grd.nz+1,1);               % tracer mixing coefficient (m2 s-1)
Mmntm.small = 1e-12;                        % arbitrarily small value
Mmntm.U = zeros(Grd.nz,1);                  % U velocity (m/s)    
Mmntm.V = zeros(Grd.nz,1);                  % V velocity (m/s)
Mmntm.q2 = ones(Grd.nz+1,1)*Mmntm.small;    % turbulent energy
Mmntm.q2l = ones(Grd.nz+1,1)*Mmntm.small;   % turbulent energy x length scale
Mmntm.len = ones(Grd.nz+1,1)*Mmntm.small;   % Mellor-Yamada 2.5 parameters
Mmntm.gh = zeros(Grd.nz+1,1);               % MY 2.5 intermediate quantity
Mmntm.sh = zeros(Grd.nz+1,1);               % MY 2.5 intermediate quantity 
Mmntm.sm = zeros(Grd.nz+1,1);               % MY 2.5 intermediate quantity
Mmntm.Km = zeros(Grd.nz+1,1);               % Momentum mixing (m2 s-1)
Mmntm.Kq = zeros(Grd.nz+1,1);               % mixing of turbulent quantities (m2 s-1)
Mmntm.boygr = zeros(Grd.nz+1,1);            % buoyancy generation term
Mmntm.shear = zeros(Grd.nz+1,1);            % shear term
Mmntm.lc_q2 = zeros(Grd.nz+1,1);            % MY 2.5 param
Mmntm.lc_q2l = zeros(Grd.nz+1,1);           % MY 2.5 param
Mmntm.kmol = In.kmol;

%--------------------------
% Calculate surface wind 
% stress parameters   
%--------------------------

% calculate wind speed at 10m (Wnd.wspeed10), wind stress (wtau) and the
% shear velocity (wstar).  The quantities are all used in approximating
% mixing coefficients and for the heat flux calculations.  The wstress
% routine comes from Rich Signell's RPSstuff toolbox.  This gives the
% stress using the Large and Pond (1981) formula and calculates the 10m 
% wind velocity based on an assumed logarithmic wind velocity profile in a
% neutral atmosphere.  This toolbox can be found at:
%
% http://woodshole.er.usgs.gov/staffpages/rsignell/rsignell.html

[ustau, vstau, uwspeed10, vwspeed10] = wstress(WndSpd.data(:,1), WndSpd.data(:,2), In.whgt);
% [Wnd.ustau,Wnd.vstau,uwspeed10,vwspeed10] = wstress(uwspeed,vwspeed,In.whgt);

wspeed10 = abs(uwspeed10 + sqrt(-1)*vwspeed10);

% Translate from dynes/cm2 to Newtons/m2

ustau = 0.1*ustau;
vstau = 0.1*vstau;

% Combine for quicker interpolation later

Wnd.t = WndSpd.t;
Wnd.data = [ustau vstau wspeed10];

%--------------------------
% Impose a pressure 
% gradient acceleration 
% for calculating mean 
% currents  
%--------------------------

% This could be made time dependent, but is held constant for now. A mean
% current, or a tidal current is crucial for ensuring some mixing in the
% interior of highly stratified systems where wind-driven currents fail to
% reach the bottom boundary

Mmntm.cor = 2*7.29e-5*sin((In.Lat/360)*2*pi); % Coriolis  

Mmntm.pgx = In.pgx;  	     % acceleration due to a pressure gradient in the EW direction
Mmntm.pgy = In.pgy;        % accelaration due to a pressure gradient in the NS direction

% Calculate bottom friction coefficient (This is the POM parameterization)
Mmntm.kappa = 0.4;
Mmntm.z0b = 0.01;
% Mmntm.Cbot = max(Mmntm.kappa^2/(log(0.5*In.dz*1/Mmntm.z0b))^2,0.0025);
% TODO: check that replacing dz with dz of bottom layer works here
Mmntm.Cbot = max(Mmntm.kappa^2/(log(0.5*In.dz(end)*1/Mmntm.z0b))^2,0.0025);

%--------------------------
% Some calculations 
% necessary for the heat 
% fluxes                
%--------------------------

% Construct a time series of the daily mean irradiance.  This quantity will
% be used in conjunction with Ht.Qo to calculate a cloud correction factor if 
% hswitch = 2.  It is also used in the biological calculations produce PAR
% averaged over 14 daylight hours.

% If you are deriving Qs, Ql, and Qlw from other meteorological inputs an
% estimate of the clear sky irradiance will be needed.  This is based on
% the "Smithsonian Formula" from Seckel and Beaudry, as reported in Reed, 
% 1977, JPO, 7, pp. 482-485.  It is good for latitudes between 20S and 60N.

Ht.Qo = clearsky(Grd.start_date,Grd.time,In.Lat)';

% calculate the mean observed daily irradiance (watts/m2).  This will
% be used in conjunction with Ht.Qo to calculate a cloud correction
% factor if hswitch = 2.  It is also used in the biological
% calculations.

winsz = 2 * round(12*3600/In.dt);
window = ones(1, winsz)/winsz;

qi = interp1(Ht.t, Ht.data(:,1), Grd.time);
% Ht.meanQi = filter(window, 1, qi);
Ht.meanQi = smooth(Grd.time, qi, winsz, 'moving'); % better edges, and no phase shift like filter

%--------------------------
% Set up archiving
%--------------------------

% TODO: Three archive variables: beginarchive, tarch, and endarchive
% In parseinput, these need to be set up
% tarch = -1: monthly
% begin/endarchive = NaN, no cutoff

nout = length(In.tarch);

% Calculate edges of bins for each output file

binedge = cell(nout,1);
for io = 1:nout
    
    % Edges, including entire simulation
    
    if In.tarch(io) == -1
        yr = Grd.start_date(1):Grd.end_date(1);
        nyr = length(yr);
        
        dv = [kron(yr', ones(12,1)) repmat((1:12)',nyr,1) ones(nyr*12,1), zeros(nyr*12,3)];
        dn = datenum(dv);
        isout = dn < dnstart | dn > dnend;
        dn = dn(~isout);

        if dn(1) ~= dnstart % mid-month start
            dn = [dnstart; dn];
        end
        if dn(end) ~= dnend
            dn = [dn; dnend];
        end
        binedge{io} = (dn - dnstart) * 86400;
    else
        binedge{io} = 0:In.tarch(io):max(Grd.time);
        if max(binedge{io}) < max(Grd.time)
            binedge{io} = [binedge{io} max(Grd.time)];
        end
        binedge{io} = binedge{io}(:);
    end
    
    % Adjust for begin/endarchive
    
    if ~isnan(In.beginarchive(io))
        t1 = (In.beginarchive(io) - dnstart)*86400;
        isbefore = binedge{io} < t1;
        binedge{io} = [0; binedge{io}(~isbefore)];  
    end
    if ~isnan(In.endarchive(io))
        t1 = (In.endarchive(io) - dnstart)*86400;
        isafter = binedge{io} > t1;
        binedge{io} = [binedge{io}(~isafter); Grd.tmax];
    end

end

% Eliminate any requested outputs that are outside the simulation timespan

isout = cellfun(@(x) length(x)<2, binedge);
if any(isout)
    binedge = binedge(~isout);
    In.outputextension = In.outputextension(~isout);
    nout = length(binedge);
    if nout == 0
        error('None of the specified output file time ranges will result in saved data');
    end
end

% Place timesteps in archiving bins

[n,bin] = deal(cell(nout,1));
for io = 1:nout
    
    [n{io}, bin{io}] = histc(Grd.time, binedge{io});
    
    % Adjust to remove final bin (where timestep equals last bin edge)
    
    n{io}(end-1) = n{io}(end-1) + n{io}(end);
    n{io}(end) = 0;
    bin{io}(bin{io} == length(binedge{io})) = length(binedge{io}) - 1;
end

% Calculate archiving variables

for io = 1:nout
    
    % Dates
    
    Arch(io).dateedge = dnstart + binedge{io}./86400;
    Arch(io).middate = (Arch(io).dateedge(1:end-1) + Arch(io).dateedge(2:end))./2;
    
    % For averaging calcs
    
    Arch(io).fraction = 1./n{io}(bin{io});
    Arch(io).islast = [logical(diff(bin{io})) true];
    Arch(io).bin = bin{io};
    Arch(io).nbin = max(Arch(io).bin);
    
	% File names: updated for new folder-based archiving, so .nc extension is 
	% removed.
    
    if nout > 1
        Arch(io).outfile = regexprep(In.outputfile, '.nc', ['_' In.outputextension{io} '.nc']);
    else
        Arch(io).outfile = In.outputfile;
    end
    [pth, fl, ex] = fileparts(Arch(io).outfile);
    Arch(io).outfile = fullfile(pth, fl);
    Arch(io).avg = [];
    Arch(io).ncid = [];
    Arch(io).vid = [];
    
    Arch(io).iens = In.iens;

end




