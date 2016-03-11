function varargout = mixed_layer(varargin)
%MIXED_LAYER Runs a 1-D physical model of a water column
%
% mixed_layer(outputfile)
% mixed_layer(outputfile, param1, val1, param2, val2, ...)
% Input = mixed_layer(...)
% Stop = mixed_layer(..., 'stopafterinit', true)
%
% This program simulates the seasonal evolution of a 1D water column.  It
% is forced by observed solar radiation and wind forcing. Numerically, a 1D
% diffusion equation is solved implicitly, allowing for long time steps
% without loss of stability.
%
% The supporting functions for this model reside in the private directory.
% Biology-related functions are store in the biomodules and biomodules_subs
% directories; both should be added to your path.  
%
% In addition, the model relies on a few 3rd-party toolboxes, which also
% need to be added to your path:
%
%   Rich Signell's RPSstuff toolbox (wstress.m) 
%   http://woodshole.er.usgs.gov/operations/sea-mat/RPSstuff-html/index.html 
%   (Because this is a large toolbox and I'm only using one function from
%   it, the most recent version of mixed_layer includes a copy of wstress.m
%   in the private directory, so you no longer need to add this toolbox).
%   
%   Phil Morgan's seawater toolbox (sw_dens0.m and sw_smow.m)
%   http://www.marine.csiro.au/datacentre/processing.htm.  
%
%   mergestruct.m: a small utility to combine structure arrays, which I
%   renamed from its original catstruct due to filename clashes with an
%   unrelated utility
%   http://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct
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
%   REQUIRED:
%   ---------
%
%   outputfile: string, base name of output folder for simulation(s).  The
%               folder name may be modified by the outputextension option
%               (see archiving options, below).  The output folder will
%               contain, at minimum, a dimensions.nc file with temporal and
%               spatial grid details and sim0001.nc file with all other
%               output variables.  If ensemble options are used (see
%               archiving options, below, as well as runmixedlayer.m), more
%               simXXXX.nc files may be added, or replaced with
%               post-processed variable files.
%
%   MODEL GRID:
%   -----------
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
%   eyear:      ending year for simulation (will end on Dec 31 of this
%               year), or 1 x 6 date vector of ending date [1976]
%
%   PHYSICAL PARAMETERS:
%   --------------------
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
%   pgx:        acceleration due to a pressure gradient in the east-west
%               direction (m s^-2) [1e-5]
%
%   pgy:        acceleration due to a pressure gradient in the north-south
%               direction (m s^-2) [0]
%
%   kmol:       molecular diffusivity (m^2/s) [1e-4]
%
%   velocity_
%   dissipation:dissipation constant (s^-1). This term removes energy from
%               past storm events over a specified time-scale as though
%               energy was being transferred to more quiescent surrounding
%               waters. [3.858e-6, i.e. 1/(3 day)]
%
%   EXTERNAL FORCING:
%   -----------------
%
%   wind_input: Wind forcing data.  This is an n x 8 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               east-west wind speed (i.e. u), north-south wind speed (i.e.
%               v). Speeds in m/s.  By default, a dateset representing 1976
%               observations on the Scotia Shelf is used.
%
%   heat_input: Heat forcing data.  This is an n x 9 matrix, with columns
%               representing year, month, day, hour, minute, second,
%               incident solar radiation (Qi), air temperature, and dew
%               point temperature.  Radiation is in W/m^2 and all
%               temperatures are in deg C. By default, a dateset
%               representing 1976 observations on the Scotia Shelf is
%               used.
%
%   ts_input:   Initial temperature and salinity profiles for simulation.
%               Data is an n x 3 matrix with columns representing depth (m,
%               negative down), temperature (deg C), and salinity (psu).
%               By default, a dateset representing 1976 observations on the
%               Scotia Shelf are used.
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
%   srelaxtime: Timescale for salinity relaxation (s) [2592000, i.e. 30 d]
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
%               salinity, and all mixed biological state variables.  NOTE:
%               This hasn't really been tested, and probably shouldn't be
%               used unless I put some more work into it.
%
%   ARCHIVING:
%   ----------
%
%   tarch:      the archiving interval (seconds). Data is averaged over
%               tarch seconds.  Can also be one of these indicators, which
%               represent non-even archiving intervals:
%               -1: monthly
%               If tarch isn't scalar, then multiple output folders will be
%               created; beginarchive and endarchive must be the same size
%               as tarch. [86400]
%
%   beginarchive:  datenumber indicating what model date to begin recording
%               to an output file, or NaN to indicate that archiving begins
%               immediately. [NaN]
%
%   endarchive: datenumber indicating what model date to stop recording
%               output to an output fle, or NaN to indicate archiving until
%               the simulation ends. [NaN]
%
%   outputextension:    cell array of strings, same size as tarch,
%               beginarchive, and endarchive.  This string is appended to
%               the outputfile string if multiple output folders are
%               indicated by the other archiving variables.  If empty and
%               tarch is non-scalar, then numeric suffixes will be used
%               (i.e. outputfile1, outputfile2, etc.)
%
%   nens:       number of ensemble members in a set of mixed_layer runs, to
%               be saved to the same output folder.  All ensemble member
%               runs *must* use the same spatial and temporal grid,
%               archiving options, and biological module; all other
%               parameters can be varied between simulations. [1]
%
%   iens:       index of the current ensemble member.  Determines the
%               number assigned to the output file name, and concatenation
%               order if postprocessing is used (see runmixedlayer.m) [1]
%
%   stopafterinit: logical scalar.  If true, simulation is terminated after
%               the initialization process, and no forward integration is
%               performed.  This is for debugging purposes.  Using this
%               option also leads to different output if mixed_layer is
%               called with an output variable (see below) [false]
%
%   BIOLOGY:
%   --------
%
%   biofun:     Function handle to biological module.  If empty, the model
%               will run without any biology. [empty] 
%
%   *var*relax: Relaxation data for any biological state variable, where
%               *var* corresponds to the short name of the variable.
%               Format is the same as for srelax.  If not included or
%               empty, no relaxation will be done.
%
%   *var*flux:  Additional flux into (or out of) any biological state
%               variable, where *var* corresponds to the short name of the
%               variable.  Format is the same as for srelax, in units of
%               stateVariableUnit/s.
%
%   brelaxtime: Timescale for relaxation of biological state variables. (s)
%               [2592000]
%
%   openbottom: Logical scalar. The value changes the way the biological
%               tracer variables interact in the bottom cell.  If false, a
%               no-flux condition is set at the bottom; mixing and vertical
%               movement will be conservative for biological tracers.  If
%               true, the values of biological variables are held constant
%               in the bottom cell, and material sinks through the
%               bottom boundary; mass is not conserved under these
%               conditions. [false]
%
%   OTHER
%   -----
%
%   verbose:    logical scalar.  If true, progress statements will be
%               printed to the screen.  If false, nothing will be printed.
%               [true]
%
%   RESTART
%   -------
%
%   Note: these options are a bit kludgy; hot starts should only be run
%   when just changing the forcing, but keeping most input (specifically,
%   archiving variables) the same.
%
%   hotstartdn: scalar datenumber, indicating time point at which
%               conditions should be saved, for use in restarting a
%               simulation at that point []
%
%   hotstartsave: name of .mat file where hot start data will be saved
%               (.mat extension should be included).  Only valid when a
%               value is supplied for hotstartdn. ['mlhotstart.mat']
%
%   hotstartload: name of file to use when restarting.  Simulation will
%               start at the last time index in the file.  All output will
%               be the same as it would have been had the simulation been
%               run from the beginning.
%
%   DEPRECATED (keeping here as a reminder, but don't use these)
%   ----------
%
%   tempfilesz: number of data points (i.e. time steps) to read in at a
%               time when converting from the temporary binary output file
%               to a netcdf file.  If NaN, all data is read at once.  Using
%               smaller amounts may speed up the read/write process, and
%               can avoid memory errors, especially if a simulation returns
%               a large number of biology-associated variables.  For
%               variable-heavy runs, reading a few hundred time steps at a
%               time is probably a good idea. [NaN]
%
%   tempdir:    string, folder where temporary file will be stored.  If
%               empty, the default temporary directory will be used. [] 
%
%   cleanup:    logical scalar.  If true, temporary files will be deleted
%               after the netcdf file is created.  Files will always be
%               kept if the simulation crashes. [true] 
%
%   newfile:    logical scalar, indicating whether to create a new set of
%               output files.  If false (as in ensemble runs), it assumes
%               the files have already been created by a previous run.
%               [true].
%
% Output variables:
%
%   Input:      1 x 1 structure holding the input variables used for the
%               run.  Returning this variable allows you to see all inputs
%               used, including those set by defaults.
%
%   Stop:       1 x 1 structure holding the variables used within a
%               mixed_layer run.  This is returned only if the
%               'stopafterinit' flag is set to true, and is useful for
%               debugging purposes.

% Copyright 2011-2015 Kelly Kearney, Charlie Stock
% kakearney@gmail.com, charles.stock@noaa.gov

%--------------------------
% Setup
%--------------------------

% Check input

In = parseinput(varargin{:});

if isfield(In, 'stopafterinit') % flag to stop after initialization
    stopafterinit = true;
    In = rmfield(In, 'stopafterinit');
else
    stopafterinit = false;
end

if In.verbose
    fprintf('\n--------------------\n1D Mixed Layer Model\n--------------------\n\n');
    fprintf('%s: %d\n\n', In.outputfile, In.iens);
end

%--------------------------
% Initialization
%--------------------------

if In.verbose
    fprintf('Initializing model...\n');
end
    
% Set constant parameters and initial conditions for simulation

[Grd, Ht, Ts, Mmntm, Wnd, Arch] = initialize(In);

% If biology is included, set initial conditions and extra parameters for
% biological calculations

if In.hasbio;
    [Bio.bio, Bio.ismixed, Bio.bottomval, Bio.Params, Bio.names, ...
        Bio.diagnames] = feval(In.biofun, 'init', In, Grd);
    Bio.ndiag = size(Bio.diagnames,1);
    if Bio.ndiag > 0
        Bio.diag = zeros(Grd.nz,Bio.ndiag); % placeholder KAK note: Charlie's new version runs biofun here to get ititial conditions
    end
    
%     Con = bioconservecheck(In.outputfile);
%     Con.dz = -diff(Grd.zp);
end

% Set up relaxation of biological variables (and additional outside fluxes)

if In.hasbio
    Bio = initbiorelax(Grd, Bio, In);
end


% Stop-after-init (for troubleshooting purposes)

if stopafterinit
    stoptempfile = tempname;
    save(stoptempfile);         % Easiest way to gather all variables into 
    A = load(stoptempfile);     % a structure is to save and reload
    delete([stoptempfile '.mat']);
    varargout{1} = A;
    if In.verbose
       fprintf('Initialization-only run; exiting mixed_layer\n');
    end
    return
end

% Set up print strings for screen display

simdatestr = datestr(datenum(Grd.start_date) + Grd.time./86400, '  mm/dd/yyyy\n');
dummystr = [repmat(' ', 1, size(simdatestr,2))-1 '\n'];
erasestr = repmat('\b', 1, size(simdatestr,2));

%--------------------------
% Mixing simulation
%--------------------------

mltimer = tic;

if ~isempty(In.hotstartload)
    Tmp = load(In.hotstartload);
    
    for ia = 1:length(Arch)
        Arch(ia).avg = Tmp.Arch(ia).avg;
    end
    
    Ts.T = Tmp.Ts.T;
    Ts.S = Tmp.Ts.S;
    Ts.Sig = Tmp.Ts.Sig;
    
    Mmntm = Tmp.Mmntm;
    
    if In.hasbio
        Bio.bio = Tmp.Bio.bio;
        Bio.diag = Tmp.Bio.diag;
    end
    
    tidx = (Tmp.Grd.savehot+1):Grd.nt;
else
    tidx = 1:Grd.nt;
end    

for it=tidx
    
    %--------------------------
    % Interpolate data to 
    % current time step
    %--------------------------
    
    tmp = interp1(Ht.t, Ht.data, Grd.time(it));
    Ht.Qi     = tmp(1);
    Ht.airtmp = tmp(2);
    Ht.dewpt  = tmp(3);
    
    tmp = interp1(Wnd.t, Wnd.data, Grd.time(it));
    Wnd.ustau    = tmp(1);
    Wnd.vstau    = tmp(2);
    Wnd.wspeed10 = tmp(3);
    
    if In.hassrelax
        Ts.srelax = interp2(Ts.Srelax.t, Ts.Srelax.z, Ts.Srelax.data, Grd.time(it), Grd.z);
    end
        
    if In.hastrelax
        Ts.trelax = interp2(Ts.Trelax.t, Ts.Trelax.z, Ts.Trelax.data, Grd.time(it), Grd.z);
    end
    
    if In.hasbio
        nb = size(Bio.names,1);
        Bio.brelax = cell(nb,1);
        for ib = 1:nb
            if Bio.hasrelax(ib)
                Bio.brelax{ib} = interp2(Bio.Relax(ib).t, Bio.Relax(ib).z, Bio.Relax(ib).data, Grd.time(it), Grd.z);
            end
        end
    end
    
    %--------------------------
    % Calculate some heat and
    % wind parameters
    %--------------------------

    % Calculate incident, sensible, latent, and longwave heat fluxes, based
    % on heat forcing/stratification at the start of each time step
     
    [Ht.Qi, Ht.Qs, Ht.Ql, Ht.Qlw, adv_fac] = calcheat(Ht.Qi, ...
            Ht.airtmp, Ht.dewpt, Ts.T(1), Wnd.wspeed10, ...
            Ht.Qo(it), Ht.meanQi(it), In.alb);
        
    % Incident heat flux is distributed throughout the water column. prad
    % is the fraction that is photosynthetically-active, attenuated
    % according to krad1.  The remaining fraction is attenuated according
    % to krad2.  Qi is also adjusted downward based on albedo.
    
    Ht.solhflx = -(1-In.alb) * Ht.Qi * diff(...
        In.prad1*exp(In.krad1*Grd.zp) + (1-In.prad1)*exp(In.krad2*Grd.zp));
    
    % Sensible, latent, and longwave fluxes are applied at the surface
    
    Ht.srfhflx = Ht.Qs + Ht.Ql + Ht.Qlw;
    
    % Advective heating/cooling fluxes
    
    if In.hasadvheat
        Ht.advhflx = In.advheatfun(Grd, Ts, it);
    else
        Ht.advhflx = zeros(size(Grd.z));
    end
        
    % Heat fluxes (W/m^2) are converted to temperature fluxes (deg C/s)
    
    Cp = 4125; % Typical specific heat for water (J kg^-1 degC^-1) @10 deg., 32 ppt)
    
    Ht.sol_Tflx = (Ht.solhflx + Ht.advhflx)./(Ts.Sig    .* Cp .* In.dz);
    Ht.srf_Tflx = Ht.srfhflx./(Ts.Sig(1) .* Cp .* In.dz(1));
   
    % Calculate bottom stress

    Wnd.ubtau = Ts.Sig(Grd.nz)*Mmntm.Cbot*sqrt(Mmntm.U(Grd.nz)^2 + ...
                Mmntm.V(Grd.nz)^2)*Mmntm.U(Grd.nz);
    Wnd.vbtau = Ts.Sig(Grd.nz)*Mmntm.Cbot*sqrt(Mmntm.U(Grd.nz)^2 + ...
                Mmntm.V(Grd.nz)^2)*Mmntm.V(Grd.nz);

    %--------------------------
    % Archiving process
    %--------------------------
                      
    % Archived variables description: variable, short name, long name,
    % units
    
    datatoarchive = {...
        Ts.T                        'temp'      'temperature'                               'deg C'
        Ts.S                        'sal'       'salinity'                                  'psu'
        Ts.Sig                      'sig'       'density'                                   'kg m^-3'
        Mmntm.U                     'ucurrent'	'east-west current velocity'                'm s^-2'
        Mmntm.V                     'vcurrent'	'north-south current velocity'              'm s^-2'
        Mmntm.q2                    'q2'        'twice the turbulent kinetic energy'        'm^2 s^-2'
        Mmntm.q2l                   'q2l'       'turbulent kinetic energy * length scale'	'm^3 s^-2'
        Mmntm.len                   'len'       'turbulence length scale'                   'm'
        Mmntm.Kh                    'Kh'        'tracer mixing coefficient'                 'm^2 s^-1'
        Mmntm.Km                    'Km'        'vertical kinematic viscosity'              'm^2 s^-1'
        Mmntm.Kq                    'Kq'        'turbulence mixing coefficient'             'm^2 s^-1'
        Mmntm.Kh .* Mmntm.boygr     'boy'       'Kh * boygr'                                'm^2 s^-3'
        Mmntm.lc_q2 .*Mmntm.q2      'diss_q2'	'lc_q2 * q2'                                '?'
        Mmntm.lc_q2l .*Mmntm.q2l    'diss_q2l'	'lc_q2l * q2l'                              '?'
        Mmntm.gh                    'gh'        'Richardson number'                         'no units'
        Ht.Qi                       'Qi'        'incident heat flux'                        'W m^-2'
        Ht.Qs                       'Qs'        'sensible heat flux'                        'W m^-2'
        Ht.Ql                       'Ql'        'latent heat flux'                          'W m^-2'
        Ht.Qlw                      'Qlw'       'longwave heat flux'                        'W m^-2'
        Ht.solhflx                  'Qsol'      'incident heat flux (water column)'         'W m^-2'
        Ht.srfhflx                  'Qsrf'      'surface heat flux'                         'W m^-2'
        Ht.advhflx                  'Qadv'      'advective heat flux'                       'W m^-2'
        Wnd.ustau                   'ustau'     'east-west surface wind stress'             'N m^-2'
        Wnd.vstau                   'vstau'     'north-south surface wind stress'           'N m^-2'
        Wnd.ubtau                   'ubtau'     'east-west bottom wind stress'              'N m^-2'
        Wnd.vbtau                   'vbtau'     'north-south bottom wind stress'            'N m^-2'
        };

               
    % If running biology, add state variables and diagnostics to this
    % matrix
    
    if In.hasbio
        if Bio.ndiag > 0
            biodata = num2cell([Bio.bio Bio.diag],1)';
            datatoarchive = [datatoarchive; [biodata [Bio.names; Bio.diagnames]]];
        else
            biodata = mat2cell(Bio.bio, size(Bio.bio,1), ones(1, size(Bio.bio,2)))';
            datatoarchive = [datatoarchive; [biodata Bio.names]];
        end
    end
               
    % Archiving
    
    for io = 1:length(Arch)
        [Arch(io).avg, Arch(io).ncid, Arch(io).vid] = archivemldata(Grd, Arch(io), it, datatoarchive, tidx(1));
    end
    if it == tidx(1)
        cu = onCleanup(@() closefiles([Arch.ncid]));
    end

     
    % Print simulation date to screen
    
    if it == tidx(1)
		if In.verbose
	        fprintf('Running simulation...\n');
	        fprintf(dummystr);
		end 
    end
    
    if Grd.newday(it) 
		if In.verbose
	        fprintf([erasestr simdatestr(it,:)]);
		end
    end

    %--------------------------
    % Calculate mixing, 
    % sinking, and growth over 
    % this time step
    %--------------------------
    
    % Calculate Mellor Yamada mixing parameters and coefficients
    % These are essentially empirical.
  
    [Mmntm.b1, Mmntm.sh, Mmntm.sm, Mmntm.Kq, Mmntm.Km, Mmntm.Kh] = ...
        mycoef(Grd.nz, Mmntm.sh, Mmntm.gh, Mmntm.sm, Mmntm.len, ...
               Mmntm.q2, Mmntm.kmol);
        

    % Calculate changes in turbulent quantities 

    [Mmntm.boygr, Mmntm.shear, Mmntm.len, Mmntm.gh, Mmntm.q2, Mmntm.q2l, ...
        Mmntm.lc_q2, Mmntm.lc_q2l] = mixturb_my(Wnd.ustau, Wnd.vstau, ... 
                    Wnd.ubtau, Wnd.vbtau, Ts.Sig, Grd.z, ... 
                    Mmntm.b1, Mmntm.U, Mmntm.V, In.dz, In.dt, Mmntm.gh, ...
                    Mmntm.len, Mmntm.Km, Mmntm.lc_q2, Mmntm.q2, Mmntm.Kh, ...
                    Mmntm.Kq, Grd.zp, In.zbot, Mmntm.q2l);
                    
    % Mix temperature and salinity 

    Ts.T = mixtracer(Ts.T, Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                     'sflux', Ht.srf_Tflx, 'bflux', 0, ...
                     'source', Ht.sol_Tflx);
                 
    Ts.S = mixtracer(Ts.S, Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                     'sflux', 0, 'bflux', 0);
                 
    % Relax temperature and salinity

    if In.hastrelax
        Ts.T = Ts.T + (Ts.trelax - Ts.T)*In.dt./In.trelaxtime;
    end
    
    if In.hassrelax
        Ts.S = Ts.S + (Ts.srelax - Ts.S)*In.dt./In.srelaxtime;
    end
    
    % Relax biology
    
    if In.hasbio
        for ib = 1:nb
            if Bio.hasrelax(ib)
                leavealone = isnan(Bio.brelax{ib});
                Bio.bio(~leavealone,ib) = Bio.bio(~leavealone,ib) + (Bio.brelax{ib}(~leavealone) - Bio.bio(~leavealone,ib))*In.dt/In.brelaxtime;
            end
        end
    end
    
    % Calculate density

    Ts.Sig = sw_dens0(Ts.S, Ts.T);

    % Advection and diffusion of water currents

	[Mmntm.U Mmntm.V] = solve_velocities([Mmntm.U; Mmntm.V],Mmntm.Km(1:Grd.nz), ...
	                    In.dt,In.dz,Mmntm.cor,Mmntm.pgx,Mmntm.pgy, ...
	                    'sflux_u',Wnd.ustau/(Ts.Sig(1)*In.dz(1)), ...
	                    'sflux_v',Wnd.vstau/(Ts.Sig(1)*In.dz(1)), ...
	                    'bflux_u',Wnd.ubtau/(Ts.Sig(end)*In.dz(end)), ...
	                    'bflux_v',Wnd.vbtau/(Ts.Sig(end)*In.dz(end)), ...
	                    'dissipate',In.velocity_dissipation);
                    

    % Mix biological variables

    if In.hasbio
        
        % Turbulent mixing
        
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 1);
            
        premixbio = Bio.bio;
        for ibio = 1:size(Bio.bio,2)
            if Bio.ismixed(ibio)
                if isnan(Bio.bottomval(ibio))
                    Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
                                    Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                                    'bflux', 0, 'sflux', 0);
                else
                    Bio.bio(:,ibio) = mixtracer(Bio.bio(:,ibio), ...
                                    Mmntm.Kh(1:Grd.nz), In.dt, In.dz, ...
                                    'bval', Bio.bottomval(ibio), ...
                                    'sflux', 0);
                end
            end
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 2);
            
        % NaN check
        
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
        % Open bottom indicator (I do this here rather than during initial
        % setup due to some issues that arose in my bio models when state
        % variables were 0 at the bottom.  By waiting until after one
        % mixing calculation, I ensure that all bio has at least a
        % minuscule amount of biomass.  A bit of a hack, and may not work
        % in all cases, so might need to look into this further later).

        if it == 1 && In.openbottom
            Bio.bottomval = Bio.bio(end,:);
        end
        
        % Other vertical movement
        
        PhysParams = Ts;
        PhysParams.par = Ht.Qi.*In.prad1;
        PhysParams.par24 = Ht.meanQi(it).*In.prad1;
        PhysParams.kpar = In.krad1;
        
        GrdParams = struct('z', Grd.z, 'dz', In.dz, 't', Grd.time(it), 'dt', In.dt);
                
        Bio.wsink = feval(In.biofun, 'vertmove', Bio.bio, ...
                              PhysParams, Bio.Params, GrdParams);
            
        presinkbio = Bio.bio;
        for ibio = 1:size(Bio.bio,2)
            Bio.bio(:,ibio) = verticalflux(Bio.bio(:,ibio), Bio.wsink(:,ibio), In.dt, In.dz, In.openbottom);
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 3);
     
          
    end
    
    % Upwelling (or downwelling): Not conservative, still trial version,
    % don't use for now
    
    
    if In.hasw
        
        w = interp2(Ts.Upwell.t, Ts.Upwell.z, Ts.Upwell.data, Grd.time(it), Grd.z);
        
        Ts.T = verticalflux(Ts.T, w, In.dt, In.dz);
        Ts.S = verticalflux(Ts.S, w, In.dt, In.dz);
        
        if In.hasbio
            for ibio = 1:size(Bio.bio,2)
                if Bio.ismixed(ibio)
                    Bio.bio(:,ibio) = verticalflux(Bio.bio(:,ibio), w, In.dt, In.dz);
                end
            end
        end
    end
   
    % A quick check to make sure that things aren't blowing up.
  
    kin_energy(it) = sum(Mmntm.U.^2 + Mmntm.V.^2);
    if any(isnan(Mmntm.U+Mmntm.V))
        error('NaNs in velocity field');
    end

    % CAS: a quick check for imaginary numbers
    
    if any(~isreal(Mmntm.U + Mmntm.V + Ts.S + Ts.T))
        error('Imaginary numbers in physical fields');
    end
  
    % Solve for biological sources/sinks and step forward

    if In.hasbio
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         PhysParams = Ts;
%         PhysParams.par = Ht.Qi.*In.prad1;
%         PhysParams.par24 = Ht.meanQi(it).*In.prad1;
%         PhysParams.krad1 = In.krad1;
%         
%         GrdParams = struct('z', Grd.z, 'dz', In.dz, 't', Grd.time(it), 'dt', In.dt);
%         
        [Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
                              PhysParams, Bio.Params, GrdParams);
        
%         [Bio.bio, Bio.diag] = feval(In.biofun, 'sourcesink', Bio.bio, ...
%                                   Ht.meanQi(it).*In.prad1, Ts.T, Grd.z, In.dz, ...
%                                   Bio.Params, Grd.time(it), In.dt);
                              
        if any(isnan(Bio.bio(:)))
            warning('ML: NaN in biology')
        end
        
%         Con = bioconservecheck(Con, Bio.bio, it, 4);
    end
    
    % Add additional outside fluxes
    
    if In.hasbio
        for ib = 1:nb
            if Bio.hasflux(ib)
                dbdtExtra = interp2(Bio.ExtraFlux(ib).t, Bio.ExtraFlux(ib).z, Bio.ExtraFlux(ib).data, Grd.time(it), Grd.z);
                Bio.bio(:,ib) = Bio.bio(:,ib) + dbdtExtra .* In.dt;
            end
        end
    end
    
    % Save data for future hot start
    
    if it == Grd.savehot
        save(In.hotstartsave);
    end


end

runtime = toc(mltimer);

% Print loop runtime

if In.verbose
	fprintf('\nMixed layer model completed successfully: %f s\n', runtime);
end

% Return input variables if output requested

if nargout == 1
    varargout{1} = In;
end

%--------------------------
% Cleanup
%--------------------------

% Debugging

if In.hasbio && isfield(Bio.Params, 'cfid')
    fclose(Bio.Params.cfid);
end


function closefiles(id)
for ii = 1:length(id)
    netcdf.close(id(ii));
end
