function varargout = mixed_layer(varargin)
%MLOFFLINE Decoupled way of running mixed_layer biological modules
%
% mixed_layer(outputfile)
% mixed_layer(outputfile, param1, val1, param2, val2, ...)
% Input = mixed_layer(...)
%
% TODO: work in progress, may crash
%
% This program simulates the seasonal evolution of a 1D water column.  It
% is forced by observed solar radiation and wind forcing. Numerically, a 1D
% diffusion equation is solved implicitly, allowing for long time steps
% without loss of stability.
%
% The supporting functions for this model reside in the private directory.
% In addition, the model relies on a few 3rd-party toolboxes:
%
%   Rich Signell's RPSstuff toolbox (wstress.m) 
%   http://woodshole.er.usgs.gov/operations/sea-mat/RPSstuff-html/index.html 
%   
%   Phil Morgan's seawater toolbox (sw_dens0.m and sw_smow.m)
%   http://www.marine.csiro.au/datacentre/processing.htm.  
%
%   mexnc and snctools
%   http://mexcdf.sourceforge.net/
%
% The mixed_layer model is set up to allow interchangeable biological
% modules to be run within it.  Please see biomodule.m (in the biomodules
% directory) for a template function and more information.
%
% The output filename input is required for the model to run.  All other
% inputs are optional and passed as parameter/value pairs.  One or more of
% these paramters can also be passed in structure format, where the
% fieldnames correspond to one or more of the parameter names.  The default
% values for these parameters set up the model for a specific example
% scenario.  Some variables (such as albedo and the attenuation
% coefficients) can be relied on for other models, while others (such as
% years, depths, or various input data sets) would make little sense when
% combined with other data. Unless indicated, all variables hold numerical
% arrays.  Variables are scalars unless dimensions are specified.  Default
% values are in brackets.
%
% Input variables:
%
%   outputfile: name of netcdf output file, where all model results will be
%               saved
%
%   dz:         the thickness of each grid cell (m).  Can be a vector of
%               thicknesses, prescribing the thickness of each individual
%               layer, although math is not as certain for this. [5]
%               
%   zbot:       the depth of the modeled water column (m, negative) [-150]
% 
%   dt:         the model time step (seconds) [10800]
%
%   syear:      starting year for simulation (will start on Jan 1 of this
%               year), or 1 x 6 date vector of starting date [1976]
%
%   eyear:      ending year for simulation (will end of Dec 31 of this
%               year), or 1 x 6 date vector of ending date [1976]
%
%   tarch:      the archiving interval (seconds). Data is averaged over
%               tarch seconds and saved in association with the mid-point
%               of each archiving interval. [86400]
%
%   krad1:      the attenuation coefficient (m^-1, value should be
%               positive) for visible radiation (~between 350 nm and 700 nm
%               wavelength).  This is approximately equivalent to the
%               photosynthetically available radiation (PAR). [0.15]
%
%   prad1:      the fraction of incoming solar radiation that falls into
%               visible wavelengths ~ PAR.  By default, it is assumed that
%               roughly 45% of the incoming solar radiation falls into this
%               category (Baker and Frouin, 1987, L&O, 32:6, pp.
%               1370-1377). [0.45]
%
%   krad2:      the attenuation coefficient (m^-1, value should be
%               positive) for non-visible (mainly infra-red) solar
%               radiation.  Water absorbs this radiation very quickly.
%               [1.67]
%
%   alb:        the albedo, or the fraction of incoming radiation reflected
%               from the sea surface.  The default value is 0.079, which is
%               typical for 43 N latitude. (Payne, R.E., 1972, Journal of
%               Atmospheric Sciences, 29:5, pp. 959-969). *If your heat
%               forcing was measured below the water surface, set the
%               albedo to 0. [0.079]
%
%   Lat:        latitude where simulation takes place (used for Coriolis
%               force calculations) [45]
%
%   whgt:       elevation above sea level where wind forcing data was
%               measured (m). [10]
%
%   wind_input: Wind forcing data.  This is an n x 8 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               east-west wind speed (i.e. u), north-south wind speed (i.e.
%               v). Speeds in m/s.  By default, a dateset representing 1976
%               observations on the Scotia Shelf are used.
%
%   heat_input: Heat forcing data.  This is an n x 9 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               incident solar radiation (Qi), air temperature, and dew
%               point temperature.  Radiation is in W/m^2 and all
%               temperatures are in deg C. By default, a dateset
%               representing 1976 observations on the Scotia Shelf are
%               used.
%
%  ts_input:    Initial temperature and salinity profiles for simulation.
%               Data is an n x 3 matrix with columns representing depth (m,
%               negative down), temperature (deg C), and salinity (psu).
%               By default, a dateset representing 1976 observations on the
%               Scotia Shelf are used.
%
%   biofun:     Function handle to biological module.  If empty, the model
%               will run without any biology. [empty] 
%
%   srelax:     Relaxtion data for salinity.  Columns 1-6 of this array
%               hold a year, month, day, hour, minute, and second
%               corresponding to the dates of the relaxation forcing, and
%               row 1 of the array hold depths (negative down, m).  Columns
%               1-6 in row 1 are just placeholders and will be ignored.
%               The remaining cells hold salinity data (ppt) for the given
%               times and depths towards which the modeled salinity will be
%               relaxed.  If all the relaxation data is from the same year,
%               the data will be treated as climatological data and will be
%               repeated for every simulated year.  If empty, no relaxation
%               will be done. [empty]
%
%   srelaxtime: Timescale for salinity relaxation (s) [2592000]
%
%   trelax:     Relaxation data for temperature.  Format is the same as for
%               srelax, with temperature data in deg C. If all the
%               relaxation data is from the same year, the data will be
%               treated as climatological data and will be repeated for
%               every simulated year.  If empty, no relaxation will be
%               done. [empty]
%
%   trelaxtime: Timescale for temperature relaxation (s) [2592000]
%
%   tracerw:    Vertical velocity data.  Format can be either a single
%               scalar value, indicating constant vertical advection over
%               space and time, or the same as for srelax, with velocity
%               data in m/s.  This data can be used to simulate upwelling
%               (positive) or downwelling (negative) velocities.  The
%               movement associated with this is applied to temperature,
%               salinity, and all mixed biological state variables. 
%
% Output variables:
%
%   Input:      1 x 1 structure holding the input variables used for the
%               run.  Returning this variable allows you to see all inputs
%               used, including those set by defaults.

% Copyright 2009 Kelly Kearney, Charlie Stock
% kkearney@princeton.edu, charles.stock@noaa.gov


%--------------------------
% Setup
%--------------------------

fprintf('\n--------------------\n1D Mixed Layer Model\n--------------------\n\n');

% Check input

In = parseinput(varargin{:});

fprintf('%s\n\n', In.outputfile);

% Add path to biomodules directory (assumes biomodules directory is located
% in same folder as mixed_layer.m)

userpath = path;
mlname = mfilename('fullpath');
mlpath = fileparts(mlname);
addpath(fullfile(mlpath, 'biomodules'));

%--------------------------
% Initialization
%--------------------------

fprintf('Initializing model...\n');

% Set constant parameters and initial conditions for simulation

[Grd, Ht, Ts, Mmntm, Wnd, Arch] = initialize(In);

% If biology is included, set initial conditions and extra parameters for
% biological calculations

if In.hasbio;
    [Bio.bio, Bio.ismixed, Bio.bottomval, Bio.Params, Bio.names, ...
        Bio.diagnames] = feval(In.biofun, 'init', In, Grd);
    Bio.ndiag = size(Bio.diagnames,1);
    if Bio.ndiag > 0
        Bio.diag = nan(Grd.nz,Bio.ndiag); % placeholder KAK note: Charlie's new version runs biofun here to get ititial conditions
    end
end

% Set up print strings for screen display

simdatestr = datestr(datenum(Grd.start_date) + Grd.time./86400, '  mm/dd/yyyy\n');
dummystr = [repmat(' ', 1, size(simdatestr,2))-1 '\n'];
erasestr = repmat('\b', 1, size(simdatestr,2));

%--------------------------
% Mixing simulation
%--------------------------

fprintf('Running simulation...\n');
fprintf(dummystr); 

tic;
for it=1:Grd.nt
    
    %--------------------------
    % Calculate some heat and
    % wind parameters
    %--------------------------
    
    % Calculate the heating based on forcing/stratifaction at the start of
    % the time step.  Qi, Qs, Ql, and Qlw are the incident, sensible, latent
    % and long-wave heat fluxes, respectively.
    
    Ht.Qi = Ht.heat(it,1);
    
%     [Ht.Qi, Ht.Qs, Ht.Ql, Ht.Qlw, adv_fac] = calcheat(Ht.heat(it,1), ...
%             Ht.heat(it,2), Ht.heat(it,3), Ts.T(1), Wnd.wspeed10(it), ...
%             Ht.Qo(it), Ht.meanQi(it), In.alb);
%     
%     % Change in incident irradiance through each depth bin
%     
%     Cp = 4125; % Typical specific heat for water (10 deg., 32 ppt)
%     
%     Ht.sol_Tflx = -(1-In.alb) * Ht.Qi * diff(In.prad1 * ...
%         exp(In.krad1 * Grd.zp) + (1-In.prad1) * exp(In.krad2*Grd.zp))./...
%         (Ts.Sig*Cp*In.dz(1));
%     
%     % Convert heatflux to a heating rate (deg. C per second) for 1 meter of
%     % water column.
%     
%     Ht.srf_Tflx = (Ht.Qs + Ht.Ql + Ht.Qlw + adv_fac) / (Ts.Sig(1) * Cp * In.dz(1));
% 
%     % Calculate bottom stress
% 
%     Wnd.ubtau = Ts.Sig(Grd.nz)*Mmntm.Cbot*sqrt(Mmntm.U(Grd.nz)^2 + ...
%                 Mmntm.V(Grd.nz)^2)*Mmntm.U(Grd.nz);
%     Wnd.vbtau = Ts.Sig(Grd.nz)*Mmntm.Cbot*sqrt(Mmntm.U(Grd.nz)^2 + ...
%                 Mmntm.V(Grd.nz)^2)*Mmntm.V(Grd.nz);

    %--------------------------
    % Archiving process
    %--------------------------
                      
    % Archived variables description: variable, short name, long name,
    % units
    
    datatoarchive = {...
        Ts.T                        'temp'      'temperature'                               'deg C'
        Ts.S                        'sal'       'salinity'                                  'psu'
        Ts.Sig                      'sig'       'density'                                   'kg m^-3'
        Ht.Qi                       'Qi'        'incident heat flux'                        'W m^-2'
        };

               
    % If running biology, add state variables and diagnostics to this
    % matrix
    
    if In.hasbio
        if Bio.ndiag > 0
            biodata = [Bio.bio Bio.diag];
            biodata = mat2cell(biodata, size(biodata,1), ones(1, size(biodata,2)))';
            datatoarchive = [datatoarchive; [biodata [Bio.names; Bio.diagnames]]];
        else
            biodata = mat2cell(Bio.bio, size(Bio.bio,1), ones(1, size(Bio.bio,2)))';
            datatoarchive = [datatoarchive; [biodata Bio.names]];
        end
    end
               
    % Create output file
    
    if it == 1
        Arch.avg = initarchivenc(datatoarchive, Grd, Arch, In.outputfile);
    end
    
    % Average data and save to file
    
    Arch.avg = archivedatanc(Arch.fraction(it), Arch.islast(it), ...
               Arch.bin(it), Arch.avg, In.outputfile, datatoarchive);
    
    % Print simulation date to screen
    
    if Arch.islast(it)
        fprintf([erasestr simdatestr(it,:)]);
    end

    %--------------------------
    % Calculate mixing, 
    % sinking, and growth over 
    % this time step
    %--------------------------
%     
%     % Calculate Mellor Yamada mixing parameters and coefficients
%     % These are essentially empirical.
%   
%     [Mmntm.b1, Mmntm.sh, Mmntm.sm, Mmntm.Kq, Mmntm.Km, Mmntm.Kh] = ...
%         mycoef(Grd.nz, Mmntm.sh, Mmntm.gh, Mmntm.sm, Mmntm.len, ...
%                Mmntm.q2, Mmntm.kmol);
% 
%     % Calculate changes in turbulent quantities 
% 
%     [Mmntm.boygr, Mmntm.shear, Mmntm.len, Mmntm.gh, Mmntm.q2, Mmntm.q2l, ...
%         Mmntm.lc_q2, Mmntm.lc_q2l] = mixturb_my(Wnd.ustau(it), Wnd.vstau(it), ... 
%                     Wnd.ubtau, Wnd.vbtau, Ts.Sig, Grd.z, ... 
%                     Mmntm.b1, Mmntm.U, Mmntm.V, In.dz, In.dt, Mmntm.gh, ...
%                     Mmntm.len, Mmntm.Km, Mmntm.lc_q2, Mmntm.q2, Mmntm.Kh, ...
%                     Mmntm.Kq, Grd.zp, In.zbot, Mmntm.q2l);
%                     
%     % Mix temperature and salinity 
% 
%     Ts.T = mixtracer(Ts.T, Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
%                      'sflux', Ht.srf_Tflx, 'bflux', 0, ...
%                      'source', Ht.sol_Tflx);
%                  
%     Ts.S = mixtracer(Ts.S, Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
%                      'sflux', 0, 'bflux', 0);
%                  
%     % Relax temperature and salinity
% 
%     if In.hastrelax
%         Ts.T = Ts.T + (Ts.trelax(:,it) - Ts.T)*In.dt./In.trelaxtime;
%     end
%     
%     if In.hassrelax
%         Ts.S = Ts.S + (Ts.srelax(:,it) - Ts.S)*In.dt./In.srelaxtime;
%     end

    Ts.T = Ts.trelax(:,it);
    Ts.S = Ts.srelax(:,it);
    
    % Calculate density

    Ts.Sig = sw_dens0(Ts.S, Ts.T);

    % Advection and diffusion of water currents

% 	[Mmntm.U Mmntm.V] = solve_velocities([Mmntm.U; Mmntm.V],Mmntm.Km(1:Grd.nz), ...
% 	                    In.dt,In.dz,Mmntm.cor,Mmntm.pgx,Mmntm.pgy, ...
% 	                    'sflux_u',Wnd.ustau(it)/(Ts.Sig(1)*In.dz(1)), ...
% 	                    'sflux_v',Wnd.vstau(it)/(Ts.Sig(1)*In.dz(1)), ...
% 	                    'bflux_u',Wnd.ubtau/(Ts.Sig(end)*In.dz(end)), ...
% 	                    'bflux_v',Wnd.vbtau/(Ts.Sig(end)*In.dz(end)), ...
% 	                    'dissipate',In.velocity_dissipation);
%                     
% 
%     % Mix biological variables
% 
%     if In.hasbio
%         
%         % Turbulent mixing
%         
%         for ibio = 1:size(Bio.bio,2)
%             if Bio.ismixed(ibio)
%                 if isnan(Bio.bottomval(ibio))
%                     Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
%                                     Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
%                                     'bflux', 0, 'sflux', 0);
%                 else
%                     Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
%                                     Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
%                                     'bval', Bio.bottomval(ibio), ...
%                                     'sflux', 0);
%                 end
%             end
%         end
%         
%         % Other vertical movement
%                 
%         Bio.wsink = feval(In.biofun, 'vertmove', Bio.bio, ...
%                           Ht.meanQi(it), Ts.T, Grd.z, In.dz, ...
%                           Bio.Params, Grd.time(it), In.dt);
%                       
%         for ibio = 1:size(Bio.bio,2)
%             Bio.bio(:,ibio) = verticalflux(Bio.bio(:,ibio), Bio.wsink(:,ibio), In.dt, In.dz);
% %             Bio.bio(:,ibio) = advectsemilag(Grd.z', Bio.bio(:,ibio)', Ts.wfun, Grd.time(it), In.dt, 'nearest')';
%         end
%         
%     end
%     
%     % Upwelling (or downwelling)
%     
%     
%     if In.hasw
%         Ts.T = advectsemilag(Grd.z, Ts.T, Ts.wfun, Grd.time(it), In.dt, 'nearest');
%         Ts.S = advectsemilag(Grd.z, Ts.S, Ts.wfun, Grd.time(it), In.dt, 'nearest');
%         if In.hasbio
%             for ibio = 1:size(Bio.bio,2)
%                 if Bio.ismixed(ibio)
%                     Bio.bio(:,ibio) = advectsemilag(Grd.z, Bio.bio(:,ibio), Ts.wfun, Grd.time(it), In.dt, 'nearest');
%                 end
%             end
%         end
%     end
%    
%     % A quick check to make sure that things aren't blowing up.
%   
%     kin_energy(it) = sum(Mmntm.U.^2 + Mmntm.V.^2);
%     if any(isnan(Mmntm.U+Mmntm.V))
%         nc_addhist(In.outputfile, 'Simulation crashed: NaNs in velocity field');
%         error('NaNs in velocity field');
%     end
% 
%     % CAS: a quick check for imaginary numbers
%     
%     if any(~isreal(Mmntm.U + Mmntm.V + Ts.S + Ts.T))
%         nc_addhist(In.outputfile, 'Simulation crashed: Imaginary numbers in physical fields');
%         error('Imaginary numbers in physical fields');
%     end
  
    % Solve for biological sources/sinks and step forward

    if it == 48 % DEBUGGING
       blah=1;
    end
        
        
    if In.hasbio
        [Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
                                  Ht.meanQi(it), Ts.T, Grd.z, In.dz, ...
                                  Bio.Params, Grd.time(it), In.dt);
    end
    
    % DEBUGGING
    
    test = all(isnan(Bio.bio(2:end,:)));
    if any(~Bio.ismixed & ~test')
        blah
    end
    
    
%     if it == 1  % DEBUGGING
%         hfig = figure;
%         hax = axes;
%         hline = plot(Bio.bio(:,2), Grd.z);
%         set(hax, 'ylim', [-200 0]);
%     else
%         set(hline, 'xdata', Bio.bio(:,2));
%         drawnow;
%     end

end

runtime = toc;

% Print loop runtime

fprintf('\nMixed layer model completed successfully: %f s\n', runtime);

% Add history attribute to file indicating successfull run

nc_addhist(In.outputfile, 'Mixed_layer simulation completed successfully');

% Return input variables if output requested

if nargout == 1
    varargout{1} = In;
end

%--------------------------
% Cleanup
%--------------------------

% Restore path to user's

path(userpath);