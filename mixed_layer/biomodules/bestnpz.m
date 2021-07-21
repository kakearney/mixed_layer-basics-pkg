function varargout = biomodule(action, varargin)
%BIOMODULE Template biological module for mixed layer model
%
% [bioinit, ismixed, bottomval, Biovars, names] = ...
%    biomodule('init', In, Grd);
% [newbio, diag] = ...
%    biomodule('sourcesink', oldbio, PhysParams, BioParams, GrdParams);
% wsink = ...
%    biomodule('vertmove', oldbio, PhysParams, BioParams, GrdParams);
%
% This is a template for a biological module for the mixed layer model.
% The module is called in three different modes during the mixed_layer
% simulation: 'init', 'sourcesink', and 'vertmove'.  
% 
% The initialize ('init') mode sets up several different variables
% associated with the biology of the model, including: 
%   - setting initial depth profiles of all biological state variables
%   - determining whether each state variable will be mixed via the
%     diffusive mixing scheme used in the model
%   - setting any deep water forcing of each state variable
%   - adding any additional variables that will be needed to run the 
%     'sourcesink' and/or 'vertmove' modes of the module. 
%   - setting the names of any dignostic variables that will be calculated
%     in the 'sourcesink' mode. 
%
% The source and sink mode ('sourcesink') calculates the change in biomass
% of each state variable over a single time step of the model due to any
% non-mixing-related process.  This is where the primary calculations will
% go for most models.
%
% The vertical movement ('vertmove') mode calculates speed of vertical
% movement for each state variable.  This can be used to introduce sinking
% or any other vertical movement unrelated to turbulent mixing.
%
% Input variables for 'init' mode:
%
%   In:         1 x 1 structure holding user-supplied info.  
%               
%               See each individual function in this folder for the
%               particular variables required by each module.
%
%   Grd:        1 x 1 structure, see mixed_layer code for fields
%
% Output variables for 'init' mode:
%
%   bioinit:    ndepth x nbsv array, initial depth profiles for all
%               biological state variables
%
%   ismixed:    1 x nbsv logical array, indicates whether each biological
%               state variable undergoes physical mixing.
%
%   bottomval:  1 x nbsv array, value used to force bottom layer during
%               mixing.  NaN indicates no forcing.
%
%   Biovars:    1 x 1 structure holding additional parameters
%
%   names:      1 x nbsv cell array of strings, names for each state
%               variable,  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
%   diagnames:  1 x ndiag cell array of strings, names for each diagnostic
%               variable.  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
% Input variables for 'sourcesink' and 'vertmove' mode:
%
%
%   oldbio:     ndepth x nbsv array, profiles of biological state variables
%               at current time step
%
%   PhysParams: 1 x 1 structure holding water column physical properties
%               relevant to biological calculations:
%
%               par:    photosynthetically active radiation (W m^-2)
%
%               par24:  photsyntheticaly active radiation, averaged over
%                       the previous 24 hours (W m^-2)  
%
%               kpar:   attenuation coefficient for PAR (m^-1)
%
%               T:      nz x 1 array of temperature
%
%               S:      nz x 1 array of salinity
%
%               Sig:    nz x 1 array of density
%
%   BioParams:  1 x 1 structure holding all biology-related parameters set
%               up in the 'init' phase.  The fields of this structure will
%               vary based on the specific biological module being used.
%
%   GrdParams:  1 x 1 structure holding parameters related to the
%               spatiotemporal model grid:
%
%               z:      ndepth x 1 array, depth grid for model (m)
%
%               dz:     depth grid cell size (m)
%
%               t:      current time (s)
%
%               dt:     time step (s)
%
% Output variables for 'sourcesink' mode:
%
%   newbio:     ndepth x nbsv array, profiles of biological state variables
%               after stepping forward by one time step
%
%   diag:       ndepth x ndiag array, values of each diagnostic variable at
%               each depth
%
% Output variables for 'vertmove' mode:
%
%   wsink:      ndepth x nbsv array, velocity of vertical movement, where
%               positive values indicate movement towards the surface (m/s)

% Copyright 2008-2015 Kelly Kearney

nin(1) = nargin(@init) + 1;
nin(2) = nargin(@sourcesink) + 1;
nin(3) = nargin(@vertmove) + 1;

nout(1) = nargout(@init);
nout(2) = nargout(@sourcesink);
nout(3) = nargout(@vertmove);

switch action
    case 'init'
        
        out = cell(1, nout(1));
				narginchk(nin(1), nin(1));       
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
				narginchk(nin(2), nin(2));        
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
				narginchk(nin(3), nin(3));        
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};
        
%**************************************************************************

function [bio, ismixed, bottomval, Biovars, names, diagnames] = init(In, Grd)

% State variables

names = {'NO3'      'nitrate'                   'mmol N m^-3' 
         'NH4'      'ammonium'                  'mmol N m^-3' 
         'PhS'      'small phytoplankton'       'mg C m^-3'
         'PhL'      'large phytoplankton'       'mg C m^-3'
         'MZS'      'small microzooplankton'    'mg C m^-3'
         'MZL'      'large microzooplankton'    'mg C m^-3'
         'Cop'      'small copepods'            'mg C m^-3'
         'NCaS'     'large copepods, on-shelf'  'mg C m^-3'
         'EupS'     'euphausiids, on-shelf'     'mg C m^-3'
         'NCaO'     'large copepods, off-shelf' 'mg C m^-3'
         'Eupo'     'euphausiids, off-shelf'    'mg C m^-3'
         'Det'      'detritus, slow-sinking'    'mg C m^-3'
         'DetF'     'detritus, fast-sinking'    'mg C m^-3'
         'Jel'      'jellyfish'                 'mg C m^-3'
         'Fe'       'iron'                      'umol Fe m^-3'
         'Ben'      'benthic infauna'           'mg C m^-2'
         'DetBen'   'detritus, benthic'         'mg C m^-2'
         'IcePhL'   'ice algae'                 'mg C m^-3'
         'IceNO3'   'ice nitrate'               'mmol N m^-3' 
         'IceNH4'   'ice ammonium'              'mmol N m^-3'
         };

% bio array is a mix of 3D variables and 2D.  Also note the crazy mix of
% units, reflecting those used in BESTNPZ.  Benthic 2D variables are stored
% in the deepest layer, and ice variables in the shallowest.

[~,kbot] = min(Grd.z);
[~,ktop] = max(Grd.z); 

bio = zeros(length(Grd.z), 20);
for ii = 1:15
    bio(:,ii) = interp1(In.bio_input3d(:,1), In.bio_input3d(:,ii+1), Grd.z);
end
bio(kbot,16:17) = In.bio_input2d(1:2);
bio(ktop,18:20) = In.bio_input2d(3:5);

% Pelagic are mixed, benthic and ice are not

ismixed = [true(1,15) false(1,5)];

% Closed system

bottomval = nan(1,20);

% No diags for now

diagnames = cell(0,3);

% All parameters can be modified on input

P = bering10kparams;
params = setdiff(fieldnames(P.Npz), {'TNU2','TNU4','AKT_BAK','TNUDG','Hout'},'stable');
p = inputParser;
p.KeepUnmatched = true;
for ii = 1:length(params)
    p.addParameter(params{ii}, P.Npz.(params{ii})(1));
end
p.addParameter('alphaPhS', 5.6); 
p.addParameter('alphaPhL', 2.2);
p.addParameter('k_sed1'  , 2.833);
p.addParameter('k_sed2'  ,-1.079);

p.parse(In);
Biovars = p.Results;

Biovars.zbot = -In.zbot;

Biovars.doyfun = @(t) doy(datetime(Grd.start_date) + seconds(t));


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)

Temp = P.T;

[~,kbot] = min(G.z);
[~,ktop] = max(G.z); 

dtdays = days(seconds(G.dt));

% First, set up the var indices, so I don't have to switch around
% between the 3 different sets used for pelagic, benthic, and ice
% variables in the input arrays, or worry about which options
% (benthic, ice, jelly, etc) are turned on.

iiNO3    = 1;
iiNH4    = 2;
iiPhS    = 3;
iiPhL    = 4;
iiMZS    = 5;
iiMZL    = 6;
iiCop    = 7;
iiNCaS   = 8;
iiEupS   = 9;
iiNCaO   = 10;
iiEupO   = 11;
iiDet    = 12;
iiDetF   = 13;
iiJel    = 14;
iiFe     = 15;
iiBen    = 16;
iiDetBen = 17;
iiIcePhL = 18;
iiIceNO3 = 19;
iiIceNH4 = 20;

% All state variables will be saved in two different versions:
% mass m^-3 (Bio3d) and mass m^-2 (Bio2d).  This redundancy
% makes for clearer (for human readers) code.

nz = length(G.z);

Bio3d = zeros(nz,20);
Bio2d = zeros(nz,20);

% Pelagic variables: These are originally stored in per-volume
% concentrations in each water column layer. NO3 and NH4 are in
% mmol N m^-3, Fe is in umol Fe m^-3, and the rest are in mg C
% m^-3.

Bio3d(:,iiNO3 ) = oldbio(:,iiNO3);
Bio3d(:,iiNH4 ) = oldbio(:,iiNH4);
Bio3d(:,iiPhS ) = oldbio(:,iiPhS);
Bio3d(:,iiPhL ) = oldbio(:,iiPhL);
Bio3d(:,iiMZS ) = oldbio(:,iiMZS);
Bio3d(:,iiMZL ) = oldbio(:,iiMZL);
Bio3d(:,iiCop ) = oldbio(:,iiCop);
Bio3d(:,iiNCaS) = oldbio(:,iiNCaS);
Bio3d(:,iiEupS) = oldbio(:,iiEupS);
Bio3d(:,iiNCaO) = oldbio(:,iiNCaO);
Bio3d(:,iiEupO) = oldbio(:,iiEupO);
Bio3d(:,iiDet ) = oldbio(:,iiDet);
Bio3d(:,iiDetF) = oldbio(:,iiDetF);
Bio3d(:,iiJel ) = oldbio(:,iiJel);
Bio3d(:,iiFe  ) = oldbio(:,iiFe);

for itrc = iiNO3:iiFe
  Bio2d(:,itrc) = Bio3d(:,itrc)*G.dz;
end

% Benthic variables: These are originally stored in per-area
% concentrations in each benthic layer, in mg C m^-2.  The
% benthic layers are of unknown thickness.
%
% At the moment, BESTNPZ hard-codes the number of benthic layers
% (NBL) to 1; for bookkeeping purposes, we're going to store
% benthic biomass in the bottom layer of our Bio3d/2d arrays.  If
% we ever change the number of benthic layers, we may need to
% rethink this schema.
%
% The 3d equivalent here is just for bookkeeping consistency (and 
% it's never used) so I'll just use the thickness of the bottom 
% layer for this conversion.

NBL = 1;

Bio2d(1,iiBen   ) = oldbio(1,iiBen);
Bio2d(1,iiDetBen) = oldbio(1,iiDetBen);

for itrc=iiBen:iiDetBen
  Bio3d(:,itrc) = Bio2d(:,itrc)./G.dz;
end

% TODO: Skipping ice, since mixed_layer doesn't do ice

% Initialize the rate of change, dB/dt, to 0 for all tracers.
% Same for all intermediate flux arrays.  Note that these
% fluxes will hold the 2D equivalent of all the fluxes (i.e.
% per area, rather than per volume); this makes it easier to
% keep track of things that are moving between different-sized
% layers (e.g. ice to surface layer, or benthos to water
% column)

DBio = zeros(nz,1); % Initializes entire array to 0

Gpp_NO3_PhS    = zeros(nz,1);
Gpp_NO3_PhL    = zeros(nz,1);
Gpp_NH4_PhS    = zeros(nz,1);
Gpp_NH4_PhL    = zeros(nz,1);
Gra_PhS_MZL    = zeros(nz,1);
Gra_PhL_MZL    = zeros(nz,1);
Ege_MZL_Det    = zeros(nz,1);
Gra_PhS_Cop    = zeros(nz,1);
Gra_PhL_Cop    = zeros(nz,1);
Gra_MZL_Cop    = zeros(nz,1);
Gra_IPhL_Cop   = zeros(nz,1);
Ege_Cop_DetF   = zeros(nz,1);
Gra_PhS_NCaS   = zeros(nz,1);
Gra_PhL_NCaS   = zeros(nz,1);
Gra_MZL_NCaS   = zeros(nz,1);
Gra_IPhL_NCaS  = zeros(nz,1);
Ege_NCaS_DetF  = zeros(nz,1);
Gra_PhS_NCaO   = zeros(nz,1);
Gra_PhL_NCaO   = zeros(nz,1);
Gra_MZL_NCaO   = zeros(nz,1);
Gra_IPhL_NCaO  = zeros(nz,1);
Ege_NCaO_DetF  = zeros(nz,1);
Gra_PhS_EupS   = zeros(nz,1);
Gra_PhL_EupS   = zeros(nz,1);
Gra_MZL_EupS   = zeros(nz,1);
Gra_Cop_EupS   = zeros(nz,1);
Gra_IPhL_EupS  = zeros(nz,1);
Gra_Det_EupS   = zeros(nz,1);
Gra_DetF_EupS  = zeros(nz,1);
Ege_EupS_DetF  = zeros(nz,1);
Gra_PhS_EupO   = zeros(nz,1);
Gra_PhL_EupO   = zeros(nz,1);
Gra_MZL_EupO   = zeros(nz,1);
Gra_Cop_EupO   = zeros(nz,1);
Gra_IPhL_EupO  = zeros(nz,1);
Gra_Det_EupO   = zeros(nz,1);
Gra_DetF_EupO  = zeros(nz,1);
Ege_EupO_DetF  = zeros(nz,1);
Gra_Cop_Jel    = zeros(nz,1);
Gra_EupS_Jel   = zeros(nz,1);
Gra_EupO_Jel   = zeros(nz,1);
Gra_NCaS_Jel   = zeros(nz,1);
Gra_NCaO_Jel   = zeros(nz,1);
Ege_Jel_DetF   = zeros(nz,1);
Mor_PhS_Det    = zeros(nz,1);
Mor_PhL_Det    = zeros(nz,1);
Mor_MZL_Det    = zeros(nz,1);
Mor_Cop_DetF   = zeros(nz,1);
Mor_NCaS_DetF  = zeros(nz,1);
Mor_EupS_DetF  = zeros(nz,1);
Mor_NCaO_DetF  = zeros(nz,1);
Mor_EupO_DetF  = zeros(nz,1);
Mor_Jel_DetF   = zeros(nz,1);
Res_PhS_NH4    = zeros(nz,1);
Res_PhL_NH4    = zeros(nz,1);
Res_MZL_NH4    = zeros(nz,1);
Res_Cop_NH4    = zeros(nz,1);
Res_NCaS_NH4   = zeros(nz,1);
Res_NCaO_NH4   = zeros(nz,1);
Res_EupS_NH4   = zeros(nz,1);
Res_EupO_NH4   = zeros(nz,1);
Res_Jel_NH4    = zeros(nz,1);
Rem_Det_NH4    = zeros(nz,1);
Rem_DetF_NH4   = zeros(nz,1);
Nit_NH4_NO3    = zeros(nz,1);
Gra_Det_Ben    = zeros(nz,1);
Gra_DetF_Ben   = zeros(nz,1);
Gra_PhS_Ben    = zeros(nz,1);
Gra_PhL_Ben    = zeros(nz,1);
Gra_DetBen_Ben = zeros(nz,1);
Exc_Ben_NH4    = zeros(nz,1);
Exc_Ben_DetBen = zeros(nz,1);
Res_Ben_NH4    = zeros(nz,1);
Mor_Ben_DetBen = zeros(nz,1);
Rem_DetBen_NH4 = zeros(nz,1);
Gpp_INO3_IPhL  = zeros(nz,1);
Gpp_INH4_IPhL  = zeros(nz,1);
Res_IPhL_INH4  = zeros(nz,1);
Mor_IPhL_INH4  = zeros(nz,1);
Nit_INH4_INO3  = zeros(nz,1);
Twi_IPhL_PhL   = zeros(nz,1);
Twi_INO3_NO3   = zeros(nz,1);
Twi_INH4_NH4   = zeros(nz,1);
Ver_PhS_DetBen = zeros(nz,1);
Ver_PhS_Out    = zeros(nz,1);
Ver_PhL_DetBen = zeros(nz,1);
Ver_PhL_Out    = zeros(nz,1);
Ver_Det_DetBen = zeros(nz,1);
Ver_Det_Out    = zeros(nz,1);
Ver_DetF_DetBen= zeros(nz,1);
Ver_DetF_Out   = zeros(nz,1);
Ver_NCaO_DetBen= zeros(nz,1);
Ver_NCaS_DetF  = zeros(nz,1);
Ver_NCaS_DetBen =zeros(nz,1);
Ver_NCaS_DetBen =zeros(nz,1);


%==============================================================
%  Biological Source/Sink terms.
%==============================================================

yday = B.doyfun(G.t);

% Copepod diapause is determined by time of year, based on sinking/
% rising day-of-year input parameters.  Set movement
% direction flags for on- and offshore large copepods here, and
% lower respiration rates if they're in the diapause (downward)
% phase.
% NCaS = CM = mostly C. marshallae, on-shelf
% NCaO = NC = mostly Neocalanus, off-shelf

RSNC = mod(B.RiseStart, 366.0);
RENC = mod(B.RiseEnd,   366.0);
SSNC = mod(B.SinkStart, 366.0);
SENC = mod(B.SinkEnd,   366.0);

% All 0 is the shortcut for lagging the onshelf group movement
% 1 month behind the offshelf group.  This maintains back-
% compatibility with input files that don't include the newer
% Rise/SinkCM parameters, since missing parameters are set to 0 by
% default.  The one-month lag was the hard-coded behavior before I
% added separate input parameters for the NCaS (i.e. CM) group.

defaultCMdp = ((B.RiseStartCM == 0.0) & ...
               (B.RiseEndCM   == 0.0) & ...
               (B.SinkStartCM == 0.0) & ...
               (B.SinkEndCM   == 0.0));

if (defaultCMdp)

    RSCM = mod(B.RiseStart + 30, 366.0);
    RECM = mod(B.RiseEnd   + 30, 366.0);
    SSCM = mod(B.SinkStart + 30, 366.0);
    SECM = mod(B.SinkEnd   + 30, 366.0);

else

    RSCM = mod(B.RiseStartCM, 366.0);
    RECM = mod(B.RiseEndCM,   366.0);
    SSCM = mod(B.SinkStartCM, 366.0);
    SECM = mod(B.SinkEndCM,   366.0);

end

upwardNC =   ((RSNC<RENC) & ...                               
             (yday>=RSNC & yday<=RENC)) ...                 
            |  ...                                               
            ((RSNC>RENC) & ...                               
             (yday>=RSNC |  yday<=RENC));

upwardCM =   ((RSCM<RECM) &  ...                              
             (yday>=RSCM & yday<=RECM)) ...                 
            | ...                                                
            ((RSCM>RECM) & ...                               
             (yday>=RSCM |  yday<=RECM));

downwardNC = ((SSNC<SENC) &  ...                              
             (yday>=SSNC & yday<=SENC)) ...                 
            | ...                                                
            ((SSNC>SENC) & ...                               
             (yday>=SSNC |  yday<=SENC));

downwardCM = ((SSCM<SECM) & ...                               
             (yday>=SSCM & yday<=SECM)) ...                 
            | ...                                                 
            ((SSCM>SECM) & ...                               
             (yday>=SSCM |  yday<=SECM));

if (downwardNC)
    respNC = B.respNCa * 0.1;
    eNC = 0;
else
    respNC = B.respNCa;
    eNC = B.eNCa;
end

if (downwardCM)
    respCM = B.respNCa * 0.1;
    eCM = 0;
else
    respCM = respNCa;
    eCM = eNCa;
end


%------------------------------
% Phytoplankton production
%------------------------------

LightLimS = ones(nz,1).*1.0;
NOLimS    = ones(nz,1).*1.0;
NHLimS    = ones(nz,1).*1.0;
NLimS     = ones(nz,1).*1.0;
IronLimS  = ones(nz,1).*1.0;
LightLimL = ones(nz,1).*1.0;
NOLimL    = ones(nz,1).*1.0;
NHLimL    = ones(nz,1).*1.0;
IronLimL  = ones(nz,1).*1.0;
NLimL     = ones(nz,1).*1.0;
PAR       = zeros(nz,1);

% Set top-of-layer irradiance to surface irradiance, 
% converted to E/m^2/d

watts2photons = 0.394848; % W m^-2 -> E/m^2/d
I0 = P.par .* watts2photons;

% Loop over layers, starting at surface...

for k = ktop:kbot
  
    % Iron limitation
    % (Hinckley et al., 2009, Deep Sea Res. II, v56(24))

    IronLimS = min(1.0, eps + Bio3d(k,iiFe)/(B.kfePhS + Bio3d(k,iiFe))*(B.kfePhS + B.FeCritPS)/B.FeCritPS); % unitless
    IronLimL = min(1.0, eps + Bio3d(k,iiFe)/(B.kfePhL + Bio3d(k,iiFe))*(B.kfePhL + B.FeCritPL)/B.FeCritPL);

    % Nitrogen limitation
    % (after COBALT, which uses Frost & Franzen, 1992)

    NOLimS = Bio3d(k,iiNO3)./((B.k1PhS + Bio3d(k,iiNO3)) .* (1.0 + Bio3d(k,iiNH4)/B.k2PhS));
    NHLimS = Bio3d(k,iiNH4)./(B.k2PhS + Bio3d(k,iiNH4));
    NOLimL = Bio3d(k,iiNO3)./((B.k1PhL + Bio3d(k,iiNO3)) .* (1.0 + Bio3d(k,iiNH4)/B.k2PhL));
    NHLimL = Bio3d(k,iiNH4)./(B.k2PhL + Bio3d(k,iiNH4));

    fratioS = NHLimS./(NOLimS + NHLimS);
    fratioL = NHLimL./(NOLimL + NHLimL);

    % Maximum uptake rate, carbon-specific and chl-specific
    % (Frost 1987,  Mar Ecol Prog Ser, v39)

    DrateS = B.DiS .* 10.0 .^ (B.DpS .* Temp(k)); % doublings d^-1 (temp dependent doubling rate)
    DrateL = B.DiL .* 10.0 .^ (B.DpL .* Temp(k)); % doublings d^-1

    PmaxS = (2.0 .^ DrateS - 1.0 );   % mg C production (mg C biomass)^-1 d^-1
    PmaxL = (2.0 .^ DrateL - 1.0 );   % mg C production (mg C biomass)^-1 d^-1

    PmaxSs = PmaxS.*B.ccr;    % mg C (mg chl)^-1 d^-1
    PmaxLs = PmaxL.*B.ccrPhL; % mg C (mg chl)^-1 d^-1

    % chl-a in layer

    chl = Bio3d(k,iiPhS)/B.ccr + Bio3d(k,iiPhL)./B.ccrPhL; % mg chl-a m^-3

    % Attenuation coefficient, including that due to clear water, 
    % chlorophyll, and optionally other organics/sediment/etc.
    % Citations for indended parameter sets are as follows:
    %
    % Luokos et al (1997, Deep Sea Res. Part II,v44(97)), after
    % Morel (1988, J. Geophys. Res., v93(C9))
    % k_ext = 0.0384, k_chlA = 0.0518, k_chlB = 0.428, 
    % k_chlC = 0, k_shallow = 0
    %
    % Ned Cokelet, personal communication (based on BEST 
    % cruises and following the method of Morel (1988)
    % k_ext = 0.034, k_chlA = 0.1159, k_chlB = 0.2829, 
    % k_chlC = 0, k_shallow = 0
    %
    % When used in the past, k_shallow = 2.0 (citation unknown)
    %
    % Update 7/17/2018: Changed hard-coded negative exponential to parameterized 
    % power law.  k_sed1 = 2.833, k_sed2 = -1.079, k_chlC = 0.0363 based
    % on fit of bottom depth vs satellite-derived Inherent Optical Properties 
    % (SNPP VIRRS absorption due to gelbstoff and detritus @ 443nm, 
    % entire-mission composite 2012-2018)

    if (B.k_sed2 < -9990.0)
        % Lazy way to allow old sediment function without recompiling 
        % (k_sed1 = old k_shallow here) (<-9990 just to avoid any floating point 
        % issues w/ -9999 equivalence)
        katten = B.k_ext + B.k_chlA.*chl.^B.k_chlB + B.k_chlC + B.k_sed1.*exp(B.zbot.*0.05);
    else
        katten = B.k_ext + B.k_chlA.*chl.^B.k_chlB + B.k_chlC + B.k_sed1.*(B.zbot).^B.k_sed2;
    end
  
    % Calculate light at depth levels relevant for Simpson's 
    % rule integration      

    z0 = 0;
    z2 = -G.dz;
    z1 = (z0+z2)/2;

    I1 = I0 .* exp(z1 .* katten);
    I2 = I0 .* exp(z2 .* katten);

    PAR(k) = (((z0-z1)./3 .* (I0 + 4.*I1 + I2))/(z0-z2))/watts2photons; % mean over layer, W m^-2

    % Calculate average light limitation across the layer
    % This approach has been adopted in order to properly 
    % capture surface production, even when using coarse 
    % vertical resolution, where light levels may not be linear 
    % within a layer (the usual convention of using the depth 
    % midpoint to represent the entire layer makes the 
    % assumption that the vertical resolution is high enough 
    % that one can assume linearity within a layer).

    % Light limitation (Jassby & Platt, 1976, Limnol Oceanogr, 
    % v21(4))

    LightLimS0 = tanh(B.alphaPhS .* I0/PmaxSs);
    LightLimS1 = tanh(B.alphaPhS .* I1/PmaxSs);
    LightLimS2 = tanh(B.alphaPhS .* I2/PmaxSs);

    LightLimL0 = tanh(B.alphaPhL .* I0/PmaxLs);
    LightLimL1 = tanh(B.alphaPhL .* I1/PmaxLs);
    LightLimL2 = tanh(B.alphaPhL .* I2/PmaxLs);

    % Light at bottom of this layer is the top of next layer

    I0 = I2;

    % Nitrate uptake, small

    f0 = Bio3d(k,iiPhS) .* PmaxS .* min([LightLimS0, NOLimS, IronLimS]);
    f1 = Bio3d(k,iiPhS) .* PmaxS .* min([LightLimS1, NOLimS, IronLimS]);
    f2 = Bio3d(k,iiPhS) .* PmaxS .* min([LightLimS2, NOLimS, IronLimS]);   

    Gpp_NO3_PhS(k) = f1; % mg C m^-3 d^-1

    % Ammonium uptake, small

    f0 = Bio3d(k,iiPhS) .* PmaxS .* min(LightLimS0, NHLimS);
    f1 = Bio3d(k,iiPhS) .* PmaxS .* min(LightLimS1, NHLimS);
    f2 = Bio3d(k,iiPhS) .* PmaxS .* min(LightLimS2, NHLimS);

    Gpp_NH4_PhS(k) = f1; % mg C m^-3 d^-1

    % Nitrate uptake, large

    f0 = Bio3d(k,iiPhL) .* PmaxL .* min([LightLimL0, NOLimL, IronLimL]);
    f1 = Bio3d(k,iiPhL) .* PmaxL .* min([LightLimL1, NOLimL, IronLimL]);
    f2 = Bio3d(k,iiPhL) .* PmaxL .* min([LightLimL2, NOLimL, IronLimL]);    

    Gpp_NO3_PhL(k) = f1; % mg C m^-3 d^-1

    % Ammonium uptake, large

    f0 = Bio3d(k,iiPhL) .* PmaxL .* min(LightLimL0, NHLimL);
    f1 = Bio3d(k,iiPhL) .* PmaxL .* min(LightLimL1, NHLimL);
    f2 = Bio3d(k,iiPhL) .* PmaxL .* min(LightLimL2, NHLimL);

    Gpp_NH4_PhL(k) = f1; % mg C m^-3 d^-1

    % Convert intermediate fluxes from volumetric to per area

    Gpp_NO3_PhS(k) = Gpp_NO3_PhS(k) .* G.dz; % mg C m^-2 d^-1
    Gpp_NH4_PhS(k) = Gpp_NH4_PhS(k) .* G.dz; % mg C m^-2 d^-1
    Gpp_NO3_PhL(k) = Gpp_NO3_PhL(k) .* G.dz; % mg C m^-2 d^-1
    Gpp_NH4_PhL(k) = Gpp_NH4_PhL(k) .* G.dz; % mg C m^-2 d^-1

    % If doubling rate is 0, ignore the above (would be more
    % efficient to check before calculating, but then I would
    % have to split up the large/small calcs; doing it this way
    % for the sake of maintenance and debugging)

    if (B.DiS <= 0.0)
        Gpp_NO3_PhS(k) = 0; % mg C m^-2 d^-1
        Gpp_NH4_PhS(k) = 0; % mg C m^-2 d^-1
    end
    if (B.DiL <= 0.0)
        Gpp_NO3_PhL(k) = 0; % mg C m^-2 d^-1
        Gpp_NH4_PhL(k) = 0; % mg C m^-2 d^-1
    end
end

%------------------------------
% Grazing and predation
%------------------------------

[BasMetCM, BasMetCop, BasMetNC, BasMetEup] = deal(zeros(nz,1));

for k=1:nz

    IcePhlAvail = 0.0;

    % Microzooplankton

    cff1 = B.fpPhSMZL .* Bio3d(k,iiPhS).^2 + ... 
         B.fpPhLMZL .* Bio3d(k,iiPhL).^2;
    cff2 = B.eMZL .* Bio3d(k,iiMZL) ./ (B.fMZL + cff1);
    cff3 = B.Q10MZL.^((Temp(k)-B.Q10MZLT)./10.0);

    Gra_PhS_MZL(k) = B.fpPhSMZL .* (Bio3d(k,iiPhS).^2) .* cff2 .* cff3; % mg C m^-3
    Gra_PhL_MZL(k) = B.fpPhLMZL .* (Bio3d(k,iiPhL).^2) .* cff2 .* cff3;

    Ege_MZL_Det(k) = (1.0 - B.gammaMZL) .* cff1 .* cff2 .* cff3; % mg C m^-2

    % Copepods

    cff1 = B.fpPhSCop .* Bio3d(k,iiPhS).^2 ...
       + B.fpPhLCop .* Bio3d(k,iiPhL).^2 ...
       + B.fpMZLCop .* Bio3d(k,iiMZL).^2 ...
       + B.fpPhLCop .* (IcePhlAvail).^2;

    cff2 = B.eCop .* Bio3d(k,iiCop) ./ (B.fCop + cff1);
    cff3 = B.Q10Cop.^((Temp(k)-B.Q10CopT)/10.0);

    if (cff1 < 0.01)  % Starvation response, used in Res below
        BasMetCop(k) = B.respCop.*cff1./0.01;
    else
        BasMetCop(k) = B.respCop;
    end

    Gra_PhS_Cop(k)  = B.fpPhSCop .* (Bio3d(k,iiPhS).^2) .* cff2 .* cff3;
    Gra_PhL_Cop(k)  = B.fpPhLCop .* (Bio3d(k,iiPhL).^2) .* cff2 .* cff3;
    Gra_MZL_Cop(k)  = B.fpMZLCop .* (Bio3d(k,iiMZL).^2) .* cff2 .* cff3;
    Gra_IPhL_Cop(k) = B.fpPhLCop .* (IcePhlAvail).^2  .* cff2 .* cff3;

    Ege_Cop_DetF(k) = (1.0 - B.gammaCop) .* cff1 .* cff2 .* cff3;

    % On-shelf Neocalanus

    cff1 = B.fpPhSNCa .* Bio3d(k,iiPhS).^2 ...
       + B.fpPhLNCa .* Bio3d(k,iiPhL).^2 ...
       + B.fpMZLNCa .* Bio3d(k,iiMZL).^2 ...
       + B.fpPhLNCa .* (IcePhlAvail).^2;

    cff2 = B.eNCa .* Bio3d(k,iiNCaS) ./ (B.fNCa + cff1);
    cff3 = B.Q10NCa .^ ((Temp(k)-B.Q10NCaT)./10.0);

    if (cff1 < 0.01) % Starvation response, used in Res below
        BasMetNC(k) = respNC.*cff1./0.01;
        BasMetCM(k) = respCM.*cff1./0.01;
    else
        BasMetNC(k) = respNC;
        BasMetCM(k) = respCM;
    end

    Gra_PhS_NCaS(k)  = B.fpPhSNCa .* Bio3d(k,iiPhS).^2 .* cff2 .* cff3;
    Gra_PhL_NCaS(k)  = B.fpPhLNCa .* Bio3d(k,iiPhL).^2 .* cff2 .* cff3;
    Gra_MZL_NCaS(k)  = B.fpMZLNCa .* Bio3d(k,iiMZL).^2 .* cff2 .* cff3;
    Gra_IPhL_NCaS(k) = B.fpPhLNCa .* (IcePhlAvail).^2 .* cff2 .* cff3;

    Ege_NCaS_DetF(k) = (1.0 - B.gammaNCa) .* cff1 .* cff2 .* cff3;

    % Off-shelf Neocalanus

    cff2 = B.eNCa .* Bio3d(k,iiNCaO) ./ (B.fNCa + cff1);

    Gra_PhS_NCaO(k)  = B.fpPhSNCa .* Bio3d(k,iiPhS).^2 .* cff2 .* cff3;
    Gra_PhL_NCaO(k)  = B.fpPhLNCa .* Bio3d(k,iiPhL).^2 .* cff2 .* cff3;
    Gra_MZL_NCaO(k)  = B.fpMZLNCa .* Bio3d(k,iiMZL).^2 .* cff2 .* cff3;
    Gra_IPhL_NCaO(k) = B.fpPhLNCa .* (IcePhlAvail).^2 .* cff2 .* cff3;

    Ege_NCaO_DetF(k) = (1.0 - B.gammaNCa) .* cff1 .* cff2 .* cff3;

    % On-shelf euphausiids

    cff1 = B.fpPhSEup * Bio3d(k,iiPhS).^2 ...
       + B.fpPhLEup * Bio3d(k,iiPhL).^2 ...
       + B.fpMZLEup * Bio3d(k,iiMZL).^2 ...
       + B.fpCopEup * Bio3d(k,iiCop).^2 ...
       + B.fpPhLEup * (IcePhlAvail).^2;    % live food

    cff0 = B.fpDetEup * Bio3d(k,iiDet).^2 ...  
       + B.fpDetEup * Bio3d(k,iiDetF).^2; % detrital food

    cff2 = B.eEup .* Bio3d(k,iiEupS) ./ (B.fEup + cff1 + cff0);
    cff3 = B.Q10Eup .^ ((Temp(k)-B.Q10EupT) ./ 10.0);

    cff4 = 1.0;

    if (cff1 < 0.01)
        BasMetEup(k) = B.respEup*cff1/0.01; % TODO: doesn't consider detrital food?
    else
        BasMetEup(k) = B.respEup;
    end

    Gra_PhS_EupS(k)  = B.fpPhSEup .* Bio3d(k,iiPhS).^2  .* cff2 .* cff3 .* cff4;
    Gra_PhL_EupS(k)  = B.fpPhLEup .* Bio3d(k,iiPhL).^2  .* cff2 .* cff3 .* cff4;
    Gra_MZL_EupS(k)  = B.fpMZLEup .* Bio3d(k,iiMZL).^2  .* cff2 .* cff3 .* cff4;
    Gra_Cop_EupS(k)  = B.fpCopEup .* Bio3d(k,iiCop).^2  .* cff2 .* cff3 .* cff4;
    Gra_IPhL_EupS(k) = B.fpPhLEup .* (IcePhlAvail).^2   .* cff2 .* cff3 .* cff4;
    Gra_Det_EupS(k)  = B.fpDetEup .* Bio3d(k,iiDet).^2  .* cff2 .* cff3 .* cff4;
    Gra_DetF_EupS(k) = B.fpDetEup .* Bio3d(k,iiDetF).^2 .* cff2 .* cff3 .* cff4;

    Ege_EupS_DetF(k) = ((1.0 - B.gammaEup) .* cff1 + ...
                        (1.0 - 0.3)   .* cff0) * ...
                        cff2 .* cff3 .* cff4;

    % Off-shelf euphausiids

    cff0 = B.fpDetEupO .* Bio3d(k,iiDet).^2 ...
       + B.fpDetEupO .* Bio3d(k,iiDetF).^2; % detrital food

    cff2 = B.eEup .* Bio3d(k,iiEupO) ./ (B.fEup + cff1 + cff0);

    cff4 = 1.0;


    Gra_PhS_EupO(k)  = B.fpPhSEup .* Bio3d(k,iiPhS).^2  .* cff2 .* cff3 .* cff4;
    Gra_PhL_EupO(k)  = B.fpPhLEup .* Bio3d(k,iiPhL).^2  .* cff2 .* cff3 .* cff4;
    Gra_MZL_EupO(k)  = B.fpMZLEup .* Bio3d(k,iiMZL).^2  .* cff2 .* cff3 .* cff4;
    Gra_Cop_EupO(k)  = B.fpCopEup .* Bio3d(k,iiCop).^2  .* cff2 .* cff3 .* cff4;
    Gra_IPhL_EupO(k) = B.fpPhLEup .* (IcePhlAvail).^2   .* cff2 .* cff3 .* cff4;
    Gra_Det_EupO(k)  = B.fpDetEupO .* Bio3d(k,iiDet).^2  .* cff2 .* cff3 .* cff4;
    Gra_DetF_EupO(k) = B.fpDetEupO .* Bio3d(k,iiDetF).^2 .* cff2 .* cff3 .* cff4;

    Ege_EupO_DetF(k) = ((1.0 - B.gammaEup) .* cff1 + ...
                      (1.0 - 0.3)   .* cff0) .* ...
                       cff2 .* cff3 .* cff4;

    % Jellyfish

    cff1 = B.fpCopJel .* Bio3d(k,iiCop).^2 + ...
         B.fpNCaJel .* Bio3d(k,iiNCaS).^2 + ...
         B.fpNCaJel .* Bio3d(k,iiNCaO).^2 + ...
         B.fpEupJel .* Bio3d(k,iiEupS).^2 + ...
         B.fpEupJel .* Bio3d(k,iiEupO).^2;

    cff2 = B.eJel .* Bio3d(k,iiJel) ./ (B.fJel + cff1);
    cff3= B.Q10Jele .^ ((Temp(k)-B.Q10JelTe) / 10.0);

    Gra_Cop_Jel(k)  = B.fpCopJel .* Bio3d(k,iiCop).^2  .* cff2 .* cff3;
    Gra_NCaS_Jel(k) = B.fpNCaJel .* Bio3d(k,iiNCaS).^2 .* cff2 .* cff3;
    Gra_NCaO_Jel(k) = B.fpNCaJel .* Bio3d(k,iiNCaO).^2 .* cff2 .* cff3;
    Gra_EupS_Jel(k) = B.fpEupJel .* Bio3d(k,iiEupS).^2 .* cff2 .* cff3;
    Gra_EupO_Jel(k) = B.fpEupJel .* Bio3d(k,iiEupO).^2 .* cff2 .* cff3;

    % Note: mentioned in Gibson & Spitz, 2011 that gammaJel can be >1 to allow
    % for an outside food source.  However, GG's code doesn't
    % specify how the flux to detritus might change in that
    % case (as written currently, that extra would come out of
    % the DetF biomass via a negative egestion flux) TODO: Do
    % we want to allow gammaJel>1, and if so, how should we
    % handle egestion?

    Ege_Jel_DetF(k) = (1.0 - B.gammaJel) .* cff1 .* cff2 .* cff3;

end

% Convert all grazing and egestion fluxes from volumetric to
% integrated over layer

for k=1:nz
    Gra_PhS_MZL(k)  = Gra_PhS_MZL(k) .* G.dz;
    Gra_PhL_MZL(k)  = Gra_PhL_MZL(k) .* G.dz;
    Ege_MZL_Det(k)  = Ege_MZL_Det(k) .* G.dz;
    Gra_PhS_Cop(k)  = Gra_PhS_Cop(k) .* G.dz;
    Gra_PhL_Cop(k)  = Gra_PhL_Cop(k) .* G.dz;
    Gra_MZL_Cop(k)  = Gra_MZL_Cop(k) .* G.dz;
    Gra_IPhL_Cop(k)  = Gra_IPhL_Cop(k) .* G.dz;
    Ege_Cop_DetF(k)  = Ege_Cop_DetF(k) .* G.dz;
    Gra_PhS_NCaS(k)  = Gra_PhS_NCaS(k) .* G.dz;
    Gra_PhL_NCaS(k)  = Gra_PhL_NCaS(k) .* G.dz;
    Gra_MZL_NCaS(k)  = Gra_MZL_NCaS(k) .* G.dz;
    Gra_IPhL_NCaS(k)  = Gra_IPhL_NCaS(k) .* G.dz;
    Ege_NCaS_DetF(k)  = Ege_NCaS_DetF(k) .* G.dz;
    Gra_PhS_NCaO(k)  = Gra_PhS_NCaO(k) .* G.dz;
    Gra_PhL_NCaO(k)  = Gra_PhL_NCaO(k) .* G.dz;
    Gra_MZL_NCaO(k)  = Gra_MZL_NCaO(k) .* G.dz;
    Gra_IPhL_NCaO(k)  = Gra_IPhL_NCaO(k) .* G.dz;
    Ege_NCaO_DetF(k)  = Ege_NCaO_DetF(k) .* G.dz;
    Gra_PhS_EupS(k)  = Gra_PhS_EupS(k) .* G.dz;
    Gra_PhL_EupS(k)  = Gra_PhL_EupS(k) .* G.dz;
    Gra_MZL_EupS(k)  = Gra_MZL_EupS(k) .* G.dz;
    Gra_Cop_EupS(k)  = Gra_Cop_EupS(k) .* G.dz;
    Gra_IPhL_EupS(k)  = Gra_IPhL_EupS(k) .* G.dz;
    Gra_Det_EupS(k)  = Gra_Det_EupS(k) .* G.dz;
    Gra_DetF_EupS(k)  = Gra_DetF_EupS(k) .* G.dz;
    Ege_EupS_DetF(k)  = Ege_EupS_DetF(k) .* G.dz;
    Gra_PhS_EupO(k)  = Gra_PhS_EupO(k) .* G.dz;
    Gra_PhL_EupO(k)  = Gra_PhL_EupO(k) .* G.dz;
    Gra_MZL_EupO(k)  = Gra_MZL_EupO(k) .* G.dz;
    Gra_Cop_EupO(k)  = Gra_Cop_EupO(k) .* G.dz;
    Gra_IPhL_EupO(k)  = Gra_IPhL_EupO(k) .* G.dz;
    Gra_Det_EupO(k)  = Gra_Det_EupO(k) .* G.dz;
    Gra_DetF_EupO(k)  = Gra_DetF_EupO(k) .* G.dz;
    Ege_EupO_DetF(k)  = Ege_EupO_DetF(k) .* G.dz;
    Gra_Cop_Jel(k)  = Gra_Cop_Jel(k) .* G.dz;
    Gra_EupS_Jel(k)  = Gra_EupS_Jel(k) .* G.dz;
    Gra_EupO_Jel(k)  = Gra_EupO_Jel(k) .* G.dz;
    Gra_NCaS_Jel(k)  = Gra_NCaS_Jel(k) .* G.dz;
    Gra_NCaO_Jel(k)  = Gra_NCaO_Jel(k) .* G.dz;
    Ege_Jel_DetF(k)  = Ege_Jel_DetF(k) .* G.dz;
end

%------------------------------
% Mortality and senescence
%------------------------------

for k = 1:nz


    % Phytoplankton (linear senescence)
    
    Mor_PhS_Det(k) = B.mPhS .* Bio3d(k,iiPhS);
    Mor_PhL_Det(k) = B.mPhL .* Bio3d(k,iiPhL);
    
    % Microzooplankton (quadratic mortality, with option for
    % linear)
    
    Mor_MZL_Det(k) = B.mpredMZL.*Bio3d(k,iiMZL).^2;   % quadratic
    
    
    TFEup = B.Q10Eup .^ ((Temp(k)-B.Q10EupT) ./ 10.0);
    
    % Mesozooplankton (quadratic predation closure)
    Mor_Cop_DetF(k)  = TFEup.*(B.mpredCop).*Bio3d(k,iiCop).^2;
    Mor_NCaS_DetF(k) = TFEup.*(B.mpredNCa).*Bio3d(k,iiNCaS).^2;
    Mor_EupS_DetF(k) = TFEup.*(B.mpredEup).*Bio3d(k,iiEupS).^2;
    Mor_NCaO_DetF(k) = TFEup.*(B.mpredNCa).*Bio3d(k,iiNCaO).^2;
    Mor_EupO_DetF(k) = TFEup.*(B.mpredEup).*Bio3d(k,iiEupO).^2;
    
    
    % Jellyfish (quadratic predation closure)
    
    Mor_Jel_DetF(k) = B.mpredJel.*Bio3d(k,iiJel).^2;

end

% Convert mortality fluxes from volumetric to integrated over
% layer

for k=1:nz

    Mor_PhS_Det(k)   = Mor_PhS_Det(k)   .* G.dz;
    Mor_PhL_Det(k)   = Mor_PhL_Det(k)   .* G.dz;
    Mor_MZL_Det(k)   = Mor_MZL_Det(k)   .* G.dz;
    Mor_Cop_DetF(k)  = Mor_Cop_DetF(k)  .* G.dz;
    Mor_NCaS_DetF(k) = Mor_NCaS_DetF(k) .* G.dz;
    Mor_EupS_DetF(k) = Mor_EupS_DetF(k) .* G.dz;
    Mor_NCaO_DetF(k) = Mor_NCaO_DetF(k) .* G.dz;
    Mor_EupO_DetF(k) = Mor_EupO_DetF(k) .* G.dz;
    Mor_Jel_DetF(k)  = Mor_Jel_DetF(k)  .* G.dz;

end

%------------------------------
% Respiration
%------------------------------

% Phytoplankton

Res_PhS_NH4 = exp(B.KtBm_PhS .* (Temp - B.TmaxPhS)) .* B.respPhS .* Bio3d(:,iiPhS);
Res_PhL_NH4 = exp(B.KtBm_PhL .* (Temp - B.TmaxPhL)) .* B.respPhL .* Bio3d(:,iiPhL);

% Microzooplankton

Res_MZL_NH4 = exp(B.KtBm_MZL .* (Temp - B.TmaxMZL)) .* B.respMZL .* Bio3d(:,iiMZL);

% Mesozooplankton (BasMetXXX is respXXX w/ starvation response)

Res_Cop_NH4  = exp(B.ktbmC .* (Temp - B.TrefC)) .* BasMetCop .* Bio3d(:,iiCop);
Res_NCaS_NH4 = exp(B.ktbmN .* (Temp - B.TrefN)) .* BasMetCM  .* Bio3d(:,iiNCaS);
Res_NCaO_NH4 = exp(B.ktbmN .* (Temp - B.TrefN)) .* BasMetNC  .* Bio3d(:,iiNCaO);
Res_EupS_NH4 = exp(B.ktbmE .* (Temp - B.TrefE)) .* BasMetEup .* Bio3d(:,iiEupS);
Res_EupO_NH4 = exp(B.ktbmE .* (Temp - B.TrefE)) .* BasMetEup .* Bio3d(:,iiEupO);

% Jellyfish (Uye & Shimauchi, 2005, J. Plankton Res. 27 (3))

Res_Jel_NH4 = B.Q10Jelr .^ ((Temp-B.Q10JelTr)./10.0) .* B.respJel .* Bio3d(:,iiJel);

% Convert respiration fluxes from volumetric to integrated over
% layer

for k=1:nz
    Res_PhS_NH4(k)  = Res_PhS_NH4(k)  .* G.dz;
    Res_PhL_NH4(k)  = Res_PhL_NH4(k)  .* G.dz;
    Res_MZL_NH4(k)  = Res_MZL_NH4(k)  .* G.dz;
    Res_Cop_NH4(k)  = Res_Cop_NH4(k)  .* G.dz;
    Res_NCaS_NH4(k) = Res_NCaS_NH4(k) .* G.dz;
    Res_NCaO_NH4(k) = Res_NCaO_NH4(k) .* G.dz;
    Res_EupS_NH4(k) = Res_EupS_NH4(k) .* G.dz;
    Res_EupO_NH4(k) = Res_EupO_NH4(k) .* G.dz;
    Res_Jel_NH4(k)  = Res_Jel_NH4(k)  .* G.dz;
end

%------------------------------
% Nitrification and
% remineralization
%------------------------------

for k=1:nz

    % Detrital remineralization

    PON = Bio3d(k,iiDet).*B.xi;  % Particulate organic nitrogen in Det, mmol N m^-3
    Rem_Det_NH4(k) = (B.Pv0 .* exp(B.PvT.*Temp(k)) .* PON); % mmol N m^-3 d^-1

    PON = Bio3d(k,iiDetF).*B.xi;  % Particulate organic nitrogen in DetF
    Rem_DetF_NH4(k) = (B.Pv0 .* exp(B.PvT*Temp(k)) .* PON); % mmol N m^-3 d^-1

    % Nitrification

    NitrifMax = B.Nitr0 .* exp(-B.ktntr.*(Temp(k) - B.ToptNtr).^2);     % Arhonditsis 2005 temperature dependence

    ParW = PAR(k); % convert to W m^-2
    DLNitrif = (1 - max(0.0, (ParW - B.tI0)./(B.KI + ParW - B.tI0))); % Fennel light dependence
    DLNitrif = 1.0;  % No light/depth dependence (overrides previous line)

    cff1 = Bio3d(k,iiNH4)./(B.KNH4Nit + Bio3d(k,iiNH4));         % Arhonditsis saturation

    Nit_NH4_NO3(k) = NitrifMax .* Bio3d(k,iiNH4) .* DLNitrif .* cff1; %  mmol N m^-3 d^-1

end

% Convert fluxes from volumetric to integrated over layer, and
% from N to C for consistency with other fluxes

for k=1:nz

    Rem_Det_NH4(k)  = Rem_Det_NH4(k)  .* G.dz./B.xi;
    Rem_DetF_NH4(k) = Rem_DetF_NH4(k) .* G.dz./B.xi;
    Nit_NH4_NO3(k)  = Nit_NH4_NO3(k)  .* G.dz./B.xi;

end


%-----------------
%Benthic Sub Model
%-----------------

% Pelagic food accessible to benthic infauna

dw = 1.0; % assume bottom 1 m is accessible

totD  = 0.0;
totDF = 0.0;
totPS = 0.0;
totPL = 0.0;

cff2 = 0.0; % accounted-for height above bottom
for k=kbot:-1:ktop

    % Fraction of this layer contributing to benthic feeding
    
    cff1 = max(min(dw, cff2 + G.dz) - cff2, 0.0); % m
    mfromlayer(k) = cff1;
    frac1(k) = cff1./G.dz;
    
    % Food available to benthos
    
    totD  = totD  + Bio2d(k,iiDet) .*frac1(k);
    totDF = totDF + Bio2d(k,iiDetF).*frac1(k);
    totPS = totPS + Bio2d(k,iiPhS) .*frac1(k);
    totPL = totPL + Bio2d(k,iiPhL) .*frac1(k);
    
    cff2 = cff2 + G.dz;
end

% Fraction of total food coming from each layer

for k=kbot:-1:ktop
    frac2(k,1) = mfromlayer(k)*Bio3d(k,iiDet) ./totD;
    frac2(k,2) = mfromlayer(k)*Bio3d(k,iiDetF)./totDF;
    frac2(k,3) = mfromlayer(k)*Bio3d(k,iiPhS) ./totPS;
    frac2(k,4) = mfromlayer(k)*Bio3d(k,iiPhL) ./totPL;
end
if (totD <= 0)
    frac2(:,1) = 0;
end
if (totDF <= 0)
    frac2(:,2) = 0;
end
if (totPS <= 0)
    frac2(:,3) = 0;
end
if (totPL <= 0)
  frac2(:,4) = 0;
end

% Potential food available from water column

cff1 = (B.prefD .*totD ./((B.prefD .*totD )+B.LupP)).*B.prefD *totD;
cff2 = (B.prefD .*totDF./((B.prefD .*totDF)+B.LupP)).*B.prefD *totDF;
cff3 = (B.prefPS.*totPS./((B.prefPS.*totPS)+B.LupP)).*B.prefPS*totPS;
cff4 = (B.prefPL.*totPL./((B.prefPL.*totPL)+B.LupP)).*B.prefPL*totPL;

cff6 = cff1+cff2+cff3+cff4; % Total pelagic food

% Potential food available from  sea floor

totBD = Bio2d(kbot,iiDetBen);
cff5 = (B.prefD .*totBD./((B.prefD .*totBD)+B.LupD)).*B.prefD *totBD;

% Temperature mediation (for feeding and mortality)

cff0 = B.q10r.^((Temp(kbot)-B.T0benr)./10.0);

% Total uptake of each food category
% TODO: Unit mismatch in the part... cff1 is mC/m^2 and
% (cff0*cff1*Bio2d(i,1,iiBen)*Rup/(cff6+KupP)) is mgC/m^2/d
% Is this supposed to be a zero-trap?  Should be cff1/dtdays
% (i.e. highest rate that would keep losses positive?)  Of
% course, that assumes no fluxes into the layer to possibly
% balance out a seemingly too-high loss rate.

cff7  = cff0.*cff1.*Bio2d(kbot,iiBen).*B.Rup./(cff6+B.KupP); % D
cff8  = cff0.*cff2.*Bio2d(kbot,iiBen).*B.Rup./(cff6+B.KupP); % DF
cff9  = cff0.*cff3.*Bio2d(kbot,iiBen).*B.Rup./(cff6+B.KupP); % PS
cff10 = cff0.*cff4.*Bio2d(kbot,iiBen).*B.Rup./(cff6+B.KupP); % PL
cff11 = cff0.*cff5.*Bio2d(kbot,iiBen).*B.Rup./(cff5+B.KupD); % DetBen

% Distribute pelagic feeding losses to appropriate water
% column layers

for k = 1:nz

    Gra_Det_Ben(k)  = cff7  .* frac2(k,1); % mg C m^-2 d^-1
    Gra_DetF_Ben(k) = cff8  .* frac2(k,2);
    Gra_PhS_Ben(k)  = cff9  .* frac2(k,3);
    Gra_PhL_Ben(k)  = cff10 .* frac2(k,4);

end

% Benthic feeding takes place in bottom layer for bookkeeping
% purposes

Gra_DetBen_Ben(kbot) = cff11; % mg C m^-2 d^-1

% Assume all excretion occurs in the bottom layer too.  Half
% goes to NH4 and half to DetBen

Exc_Ben_DetBen(kbot) = (B.eexD * (cff7 + cff8 + cff11) + ...
                       B.eex  * (cff9 + cff10)) * 0.5;
Exc_Ben_NH4(kbot) = Exc_Ben_DetBen(kbot);

% Respiration (also takes place in bottom layer)

cff3 = cff0 * Bio2d(kbot,iiBen) * B.Rres;
cff4 = ((1 - B.eexD) * (cff7 + cff8 + cff11) + ...
        (1 - B.eex)  * (cff9 + cff10)) * B.Qres;

Res_Ben_NH4(kbot) = cff3 + cff4; % mg C m^-2 d^-1

% Mortality (linear senescence and quadratic predation closure)

Mor_Ben_DetBen(kbot) = cff0.*B.rmort  .*Bio2d(kbot,iiBen) + ...
                      cff0.*B.BenPred.*Bio2d(kbot,iiBen).^2;  % mg C m^-2 d^-1

% Benthic remineralization: assumes only the top 25% is
% available to remineralize to NH4 (in bottom layer) 
% (Kawamiya et al., 2000, J. Mar. Syst., v25(2))

PON = Bio3d(kbot,iiDetBen).*0.25.*B.xi;  % Benthic Particulate organic nitrogen
cff1 = B.Pv0.*exp(B.PvT.*Temp(kbot)).*PON;  % mmol N m^-3 d^-1

Rem_DetBen_NH4(kbot) = cff1.*G.dz./B.xi; % mg C m^-2 d^-1

% Skipping more ice calcs

%------------------------------
% Combine bio source/sinks
%------------------------------

DBio(:,iiNO3   ) = (Nit_NH4_NO3                             ...
                   +  Twi_INO3_NO3                            ...
                   -  Gpp_NO3_PhS                             ...
                   -  Gpp_NO3_PhL)*B.xi*dtdays; % NO3: mmolN m^-2

DBio(:,iiNH4   ) = (Res_PhS_NH4                             ...
                   +  Res_PhL_NH4                             ...
                   +  Res_MZL_NH4                             ...
                   +  Res_Cop_NH4                             ...
                   +  Res_NCaS_NH4                            ...
                   +  Res_NCaO_NH4                            ...
                   +  Res_EupS_NH4                            ...
                   +  Res_EupO_NH4                            ...
                   +  Res_Jel_NH4                             ...
                   +  Rem_Det_NH4                             ...
                   +  Rem_DetF_NH4                            ...
                   +  Exc_Ben_NH4                             ...
                   +  Res_Ben_NH4                             ...
                   +  Rem_DetBen_NH4                          ...
                   +  Twi_INH4_NH4                            ...
                   -  Gpp_NH4_PhS                             ...
                   -  Gpp_NH4_PhL                             ...
                   -  Nit_NH4_NO3)*B.xi*dtdays; % NH4: mmol N m^-2

DBio(:,iiPhS   ) = (Gpp_NO3_PhS                             ...
                   +  Gpp_NH4_PhS                             ...
                   -  Gra_PhS_MZL                             ...
                   -  Gra_PhS_Cop                             ...
                   -  Gra_PhS_NCaS                            ...
                   -  Gra_PhS_NCaO                            ...
                   -  Gra_PhS_EupS                            ...
                   -  Gra_PhS_EupO                            ...
                   -  Mor_PhS_Det                             ...
                   -  Res_PhS_NH4                             ...
                   -  Gra_PhS_Ben)*dtdays; % PhS: mg C m^-2

DBio(:,iiPhL   ) = (Gpp_NO3_PhL                             ...
                   +  Gpp_NH4_PhL                             ...
                   +  Twi_IPhL_PhL                            ...
                   -  Gra_PhL_MZL                             ...
                   -  Gra_PhL_Cop                             ...
                   -  Gra_PhL_NCaS                            ...
                   -  Gra_PhL_NCaO                            ...
                   -  Gra_PhL_EupS                            ...
                   -  Gra_PhL_EupO                            ...
                   -  Mor_PhL_Det                             ...
                   -  Res_PhL_NH4                             ...
                   -  Gra_PhL_Ben)*dtdays; % PhL: mg C m^-2

DBio(:,iiMZL   ) = (Gra_PhS_MZL                             ...
                   +  Gra_PhL_MZL                             ...
                   -  Ege_MZL_Det                             ...
                   -  Gra_MZL_Cop                             ...
                   -  Gra_MZL_NCaS                            ...
                   -  Gra_MZL_NCaO                            ...
                   -  Gra_MZL_EupS                            ...
                   -  Gra_MZL_EupO                            ...
                   -  Mor_MZL_Det                             ...
                   -  Res_MZL_NH4)*dtdays; % MZL: mg C m^-2

DBio(:,iiCop   ) = (Gra_PhS_Cop                             ...
                   +  Gra_PhL_Cop                             ...
                   +  Gra_MZL_Cop                             ...
                   +  Gra_IPhL_Cop                            ...
                   -  Ege_Cop_DetF                            ...
                   -  Gra_Cop_EupS                            ...
                   -  Gra_Cop_EupO                            ...
                   -  Gra_Cop_Jel                             ...
                   -  Mor_Cop_DetF                            ...
                   -  Res_Cop_NH4)*dtdays; % Cop: mg C m^-2

DBio(:,iiNCaS  ) = (Gra_PhS_NCaS                            ...
                   +  Gra_PhL_NCaS                            ...
                   +  Gra_MZL_NCaS                            ...
                   +  Gra_IPhL_NCaS                           ...
                   -  Ege_NCaS_DetF                           ...
                   -  Gra_NCaS_Jel                            ...
                   -  Mor_NCaS_DetF                           ...
                   -  Res_NCaS_NH4)*dtdays; % NCaS: mg C m^-2

DBio(:,iiEupS  ) = (Gra_PhS_EupS                            ...
                   +  Gra_PhL_EupS                            ...
                   +  Gra_MZL_EupS                            ...
                   +  Gra_Cop_EupS                            ...
                   +  Gra_IPhL_EupS                           ...
                   +  Gra_Det_EupS                            ...
                   +  Gra_DetF_EupS                           ...
                   -  Ege_EupS_DetF                           ...
                   -  Gra_EupS_Jel                            ...
                   -  Mor_EupS_DetF                           ...
                   -  Res_EupS_NH4)*dtdays; % EupS: mg C m^-2

DBio(:,iiNCaO  ) = (Gra_PhS_NCaO                            ...
                   +  Gra_PhL_NCaO                            ...
                   +  Gra_MZL_NCaO                            ...
                   +  Gra_IPhL_NCaO                           ...
                   -  Ege_NCaO_DetF                           ...
                   -  Gra_NCaO_Jel                            ...
                   -  Mor_NCaO_DetF                           ...
                   -  Res_NCaO_NH4)*dtdays; % NCaO: mg C m^-2

DBio(:,iiEupO  ) = (Gra_PhS_EupO                            ...
                    +  Gra_PhL_EupO                           ...
                    +  Gra_MZL_EupO                           ...
                    +  Gra_Cop_EupO                           ...
                    +  Gra_IPhL_EupO                          ...
                    +  Gra_Det_EupO                           ...
                    +  Gra_DetF_EupO                          ...
                    -  Ege_EupO_DetF                          ...
                    -  Gra_EupO_Jel                           ...
                    -  Mor_EupO_DetF                          ...
                    -  Res_EupO_NH4)*dtdays; % EupO: mg C m^-2

DBio(:,iiDet   )  = (Ege_MZL_Det                            ...
                    +  Mor_PhS_Det                            ...
                    +  Mor_PhL_Det                            ...
                    +  Mor_MZL_Det                            ...
                    -  Gra_Det_EupS                           ...
                    -  Gra_Det_EupO                           ...
                    -  Rem_Det_NH4                            ...
                    -  Gra_Det_Ben)*dtdays; % Det: mg C m^-2

DBio(:,iiDetF  )  = (Ege_Cop_DetF                           ...
                    +  Ege_NCaS_DetF                          ...
                    +  Ege_NCaO_DetF                          ...
                    +  Ege_EupS_DetF                          ...
                    +  Ege_EupO_DetF                          ...
                    +  Ege_Jel_DetF                           ...
                    +  Mor_Cop_DetF                           ...
                    +  Mor_NCaS_DetF                          ...
                    +  Mor_EupS_DetF                          ...
                    +  Mor_NCaO_DetF                          ...
                    +  Mor_EupO_DetF                          ...
                    +  Mor_Jel_DetF                           ...
                    -  Gra_DetF_EupS                          ...
                    -  Gra_DetF_EupO                          ...
                    -  Rem_DetF_NH4                           ...
                    -  Gra_DetF_Ben)*dtdays; % DetF: mg C m^-2

DBio(:,iiJel   )  = (Gra_Cop_Jel                            ...
                    +  Gra_EupS_Jel                           ...
                    +  Gra_EupO_Jel                           ...
                    +  Gra_NCaS_Jel                           ...
                    +  Gra_NCaO_Jel                           ...
                    -  Ege_Jel_DetF                           ...
                    -  Mor_Jel_DetF                           ...
                    -  Res_Jel_NH4)*dtdays; % Jel: mg C m^-2

DBio(:,iiFe    )  = (                                       ...
                    -  Gpp_NO3_PhS                            ...
                    -  Gpp_NO3_PhL)*B.FeC*dtdays; % Fe: umol Fe m^-2

DBio(:,iiBen   )  = (Gra_Det_Ben                            ...
                    +  Gra_DetF_Ben                           ...
                    +  Gra_PhS_Ben                            ...
                    +  Gra_PhL_Ben                            ...
                    +  Gra_DetBen_Ben                         ...
                    -  Exc_Ben_NH4                            ...
                    -  Exc_Ben_DetBen                         ...
                    -  Res_Ben_NH4                            ...
                    -  Mor_Ben_DetBen)*dtdays; % Ben: mg C m^-2

DBio(:,iiDetBen)  = (Exc_Ben_DetBen                         ...
                    +  Mor_Ben_DetBen                         ...
                    -  Gra_DetBen_Ben                         ...
                    -  Rem_DetBen_NH4)*dtdays; % DetBen: mg C m^-2

DBio(:,iiIcePhL)  = (Gpp_INO3_IPhL                          ...
                    +  Gpp_INH4_IPhL                          ...
                    -  Gra_IPhL_Cop                           ...
                    -  Gra_IPhL_NCaS                          ...
                    -  Gra_IPhL_NCaO                          ...
                    -  Gra_IPhL_EupS                          ...
                    -  Gra_IPhL_EupO                          ...
                    -  Res_IPhL_INH4                          ...
                    -  Mor_IPhL_INH4                          ...
                    -  Twi_IPhL_PhL)*dtdays; % IcePhL: mg C m^-2

DBio(:,iiIceNO3)  = (Nit_INH4_INO3                          ...
                    -  Gpp_INO3_IPhL                          ...
                    -  Twi_INO3_NO3)*B.xi*dtdays; % IceNO3: mmol N m^-2

DBio(:,iiIceNH4)  = (Res_IPhL_INH4                          ...
                    +  Mor_IPhL_INH4                          ...
                    -  Gpp_INH4_IPhL                          ...
                    -  Nit_INH4_INO3                          ...
                    -  Twi_INH4_NH4)*B.xi*dtdays; % IceNH4: mmol N m^-2


% Add DBio terms to existing biomass

Bio2d = Bio2d + DBio;

% Infauna (Ben) group can receive flux from water column layers.
% Move these additions to the bottom layer now, consistent with
% the initial setup of the Bio2d and Bio3d arrays.

Bio2d(kbot,iiBen) = sum(Bio2d(:,iiBen), 2);
Bio2d(setdiff(1:nz,kbot),iiBen) = 0.0;

% TODO: Eliminate negatives?  Hopefully processes are
% formulated to prevent any, but very fast overturning might
% result in numerical issues.  Brute force zero traps will
% eliminate conservation of mass, so I'd prefer to look into
% increasing BioIter if this is a problem

Bio2d = max(Bio2d, 0.0);

% TODO: Save this "input flux" as diagnostic

% Sync volumetric version to the updated per-area values

for k = 1:nz
    for itrc = 1:17 % Pelagic (and benthic, for bookkeeping)
        Bio3d(k,itrc) = Bio2d(k,itrc)./G.dz;
    end
    for itrc = 18:20 % Ice
        Bio3d(k,itrc) = Bio2d(k,itrc)./B.aidz;
    end
end
                
% MODIFY CODE HERE
                         
newbio = oldbio;
diag = rand(size(oldbio,1), 1);

%**************************************************************************

function wsink = vertmove(oldbio, P, B, G)

% MODIFY CODE HERE

wsink = zeros(size(oldbio));


function [dBioOut, flx] = biosink(nn,wBio,Bio,HzL,dtdays,z_wL,zlimit)

