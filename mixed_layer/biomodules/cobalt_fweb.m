function varargout = cobalt_fweb(action, varargin)
%COBALT_FWEB COBALT biological module
%
% This module contains the planktonic food web components of the Carbon,
% Ocean Biogeochemistry and Lower Trophics (COBALT) model described in:
% 
% Stock, Dunne, and John, 2014. Global-scale carbon and energy flows
% through the marine planktonic food web: an analysis with a coupled
% physical-biological model.  Progress in Oceanography, 120 1-28.
%
% 12 State variables capturing the core planktonic food web interactions
% are included:
%
% sp = small phytoplankton
% lp = large phytoplankton
% sz = small zooplankton
% mz = medium zooplankton
% lz = large zooplankton
% b = bacteria
% det = detritus
% ldon = labile dissolved organic nitrogen
% sldon = semi-labile dissolved organic nitrogen
% srdon = semi-refractory dissolved organic nitrogen
% nh4 = ammonia
% no3 = nitrate
%
% See Fig. 1 of Stock et al. (2014) for the basic structure (noting that
% diazotrophs are omitted).  Other details of the parameters and their
% values can also be found in this paper.
%
% This code is provided "as is", questions or comments can be directed to
% Charles Stock (charles.stock@noaa.gov).
%
% See biomodule.m for function syntax descriptions.  The following 
% fields must be present in the In structure (passed to mixed_layer as
% parameter/value pairs):
%
%   bio_input:  Initial conditions for biological data, array has number
%               of colums equal to the number of state varials in the
%               biological module + 1, where the first column is depth
%
% Copyright 2015 Charles Stock

nin(1) = nargin(@init) + 1;
nin(2) = nargin(@sourcesink) + 1;
nin(3) = nargin(@vertmove) + 1;

nout(1) = nargout(@init);
nout(2) = nargout(@sourcesink);
nout(3) = nargout(@vertmove);

switch action
    case 'init'
        
        out = cell(1, nout(1));
        error(nargchk(nargin,nin(1),nin(1)));        
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
        error(nargchk(nargin,nin(2),nin(2)));        
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
        error(nargchk(nargin,nin(3),nin(3)));        
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};
        
%**************************************************************************

function [bio, ismixed, bottomval, params, names, diag, diagnames] = init(In, Grd)

% 12 State variables, See Fig. 1 of Stock et al. 2014 but omit diazotrophs
bio(:,1) = interp1(In.bio_input(:,1), In.bio_input(:,2), Grd.z);
bio(:,2) = interp1(In.bio_input(:,1), In.bio_input(:,3), Grd.z);
bio(:,3) = interp1(In.bio_input(:,1), In.bio_input(:,4), Grd.z);
bio(:,4) = interp1(In.bio_input(:,1), In.bio_input(:,5), Grd.z);
bio(:,5) = interp1(In.bio_input(:,1), In.bio_input(:,6), Grd.z);
bio(:,6) = interp1(In.bio_input(:,1), In.bio_input(:,7), Grd.z);
bio(:,7) = interp1(In.bio_input(:,1), In.bio_input(:,8), Grd.z);
bio(:,8) = interp1(In.bio_input(:,1), In.bio_input(:,9), Grd.z);
bio(:,9) = interp1(In.bio_input(:,1), In.bio_input(:,10), Grd.z);
bio(:,10) = interp1(In.bio_input(:,1), In.bio_input(:,11), Grd.z);
bio(:,11) = interp1(In.bio_input(:,1), In.bio_input(:,12), Grd.z);
bio(:,12) = interp1(In.bio_input(:,1), In.bio_input(:,13), Grd.z);

% State Variable Names (column vector)
names = {'sp', 'small phytoplankton conc', 'mmoles N m-3';
         'lp', 'large phytoplankton conc', 'mmoles N m-3';
         'sz', 'small zooplankton conc', 'mmoles N m-3';
         'mz', 'medium zooplankton conc', 'mmoles N m-3';
         'lz', 'large zooplankton conc', 'mmoles N m-3';
         'b', 'bacteria conc', 'mmoles N m-3';
         'det', 'detritus conc', 'mmoles N m-3';
         'ldon','labile don conc', 'mmoles N m-3';
         'sldon','semi-labile don conc','mmoles N m-3';
         'srdon','semi-refractory don conc','mmoles N m-3';
         'nh4','ammonia conc','mmoles nh4 m-3'
         'no3','nitrate conc','mmoles no3 m-3'};

% specify if mixed or not
ismixed = true(1,12);

% Nitrate is set at the bottom boundary, note that biorlx provides options
% to relax deep nutrients at multiple grid cells (not just the bottom 
% boundary) toward specified values. 
bottomval = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN bio(end,12)];
% for use this version of bottomval to check conservation
%bottomval = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

%
%   COBALT parameter values (Default if global formulation of Stock et al.,
%   2014.  Note that all rates are given for 0 deg. C
%

sperd = 86400;                      % sec day-1
% sinking
params.wsnk = -100/sperd;           % m sec-1 (negative is down)
% refuge concentration for basal metabolic rates
params.ref_conc = 1e-4;             % mmoles m-3
% Phytoplankton
params.kappa_phyto = 0.063;        % deg C-1
params.alpha_sp = 2.0e-5*2.77e18/6.022e17; % g C g Chl-1 m2 W-1 s-1
params.alpha_lp = 1.0e-5*2.77e18/6.022e17; % g C g Chl-1 m2 W-2 s-1
params.P_C_max_sp = 1.125/sperd;   % sec-1
params.P_C_max_lp = 1.25/sperd;    % sec-1
params.thetamax_sp = 0.03;         % g Chl g C-1
params.thetamax_lp = 0.05;         % g Chl g C-1
params.zeta_sp = 0.05;             % dimensionless
params.zeta_lp = 0.05;             % dimensionless
params.bresp_sp = 0.0225/sperd;    % sec-1
params.bresp_lp = 0.025/sperd;     % sec-1
params.thetamin = 0.002;           % g Chl g C-1
params.kno3_sp = 0.5;              % mmoles NO3 m-3
params.kno3_lp = 2.5;              % mmoles NO3 m-3
params.knh4_sp = 0.1;              % mmoles NH4 m-3
params.knh4_lp = 0.5;              % mmoles NH4 m-3
params.m_agg_sp = 0.1/sperd;       % sec-1 mmole N-1 m3
params.m_agg_lp = 0.3/sperd;       % sec-1 mmole N-1 m3
params.m_vir_sp = 0.025/sperd;     % sec-1 mmole N-1 m3
params.m_vir_lp = 0/sperd;         % sec-1 mmole N-1 m3
params.exu = 0.13;                 % dimensionless
params.mld_thresh = 0.01;          % kg m-3

% Zooplankton
params.kappa_zoo = 0.063;          % deg C-1
params.i_max_sz = 1.42/sperd;      % sec-1
params.i_max_mz = 0.57/sperd;      % sec-1
params.i_max_lz = 0.23/sperd;      % sec-1
params.ki_sz = 1.25;               % mmoles N m-3
params.ki_mz = 1.25;               % mmoles N m-3
params.ki_lz = 1.25;               % mmoles N m-3
params.gge_max_sz = 0.4;           % dimensionless
params.gge_max_mz = 0.4;           % dimensionless
params.gge_max_lz = 0.4;           % dimensionless
params.bresp_sz = 0.020/sperd;     % sec-1
params.bresp_mz = 0.008/sperd;     % sec-1
params.bresp_lz = 0.0032/sperd;    % sec-1

% detritus/dissolved organic matter partitioning, 30% goes to detritus,
% where it is partitioned between sinking and dissolved according to size.
% For the dissolved phase, partitioning between labile, semi-labile and
% semi-refractory is based on a tuning to global patterns.
params.phi_det_sz = 0.0;            % dimensionless
params.phi_det_mz = 0.20;            % dimensionless
params.phi_det_lz = 0.30;            % dimensionless
params.phi_ldon_sz = 0.30*0.57;      % dimensionless
params.phi_ldon_mz = 0.10*0.57;      % dimensionless
params.phi_ldon_lz = 0.0;            % dimensionless
params.phi_sldon_sz = 0.30*0.40;     % dimensionless
params.phi_sldon_mz = 0.10*0.40;     % dimensionless
params.phi_sldon_lz = 0.0;           % dimensionless
params.phi_srdon_sz = 0.30*0.03;     % dimensionless
params.phi_srdon_mz = 0.10*0.03;     % dimensionless
params.phi_srdon_lz = 0.0;           % dimensionless
% switching shape parameters, see: Stock, Powell and Levin (2008),
% Bottom-up and top-down forcing in a simple size structured plankton
% dynamics model.  Journal of Marine Systems 74(1) 134-152 for rationale
% and references.
params.ms_sz = 2;
params.ns_sz = 2;
params.ms_mz = 2;
params.ns_mz = 2;
params.ms_lz = 2;
params.ns_lz = 2;
% innate prey availability for sp, lp, sz, mz, lz, ndet; dimension nz x 7
ipa_matrix_sz = [1 0 0 0 0 0.25 0];
ipa_matrix_mz = [0 1 1 0 0 0 0];
ipa_matrix_lz = [0 1 0 1 0 0 0];
for n = 2:size(Grd.z,1)
    ipa_matrix_sz = [ipa_matrix_sz; 1 0 0 0 0 0.25 0];
    ipa_matrix_mz = [ipa_matrix_mz; 0 1 1 0 0 0 0];
    ipa_matrix_lz = [ipa_matrix_lz; 0 1 0 1 0 0 0];
end

params.ipa_matrix_sz = ipa_matrix_sz;
params.ipa_matrix_mz = ipa_matrix_mz;
params.ipa_matrix_lz = ipa_matrix_lz;


% Higher predation closure
params.kappa_hp = 0.063;           % deg C-1
params.i_max_hp = 0.09/sperd;      % sec-1
params.ki_hp = 1.25;               % mmoles N m-3
params.ms_hp = 2;           % dimensionless
params.ns_hp = 2;           % dimensionless
ipa_matrix_hp = [0 0 0 1 1 0 0];   % dimensionless
for n = 2:size(Grd.z,1)
    ipa_matrix_hp = [ipa_matrix_hp; 0 0 0 1 1 0 0];
end
params.ipa_matrix_hp = ipa_matrix_hp;
params.coef_hp = 2;                % dimensionless
% partitioning of higher predator losses to different pools
params.phi_det_hp = 0.35;          % dimensionless
params.phi_ldon_hp = 0.0;          % dimensionless
params.phi_sldon_hp = 0.0;         % dimensionless
params.phi_srdon_hp = 0.0;         % dimensionless
params.phi_nh4_hp = 0.65;          % dimensionless

% Bacteria
params.kappa_bact = 0.063;         % deg C-1
params.mu_max_b = 1.0/sperd;       % sec-1             
params.kldon = 0.5;                % mmoles ldon m-3
params.gge_max_b = 0.4;            % dimensionless
params.bresp_b = 0.0075/sperd;     % sec-1
params.m_vir_b = 0.033/sperd;      % sec-1 mmole N-1 m3
% partitioning of virus losses amongest dissolved organics
params.phi_ldon_vir = 0.55;        % dimensionless
params.phi_sldon_vir = 0.40;       % dimensionless
params.phi_srdon_vir = 0.05;       % dimensionless

% Detritus and dissolved organic material
params.gamma_det = -params.wsnk/188;        % sec-1
params.gamma_srdon = 1.0/(18*365*sperd);    % sec-1
params.gamma_sldon = 1.0/(90*sperd);        % sec-1
params.gamma_nitrif = 1.0/(30*sperd);       % sec-1

%
% Initialize diagnostic variables, whatever is added here must be
% calculated in sourcesink
%
diag.spprod = zeros(size(Grd.z));  % small phyto prod (mg N m-3 day-1)
diag.lpprod = zeros(size(Grd.z));  % large phyto prod (mg N m-3 day-1)
diag.bprod = zeros(size(Grd.z));   % bacteria prod (mg N m-3 day-1)
diag.szprod = zeros(size(Grd.z));  % small zoo prod (mg N m-3 day-1)
diag.mzprod = zeros(size(Grd.z));  % medium zoo prod (mg N m-3 day-1)
diag.lzprod = zeros(size(Grd.z));  % large zoo prod (mg N m-3 day-1)
diag.theta_sp = zeros(size(Grd.z)); % small chl:C ratio
diag.theta_lp = zeros(size(Grd.z)); % large chl:C ratio
diag.chl = zeros(size(Grd.z));      % chlorophyll (mg chl m-3)

diagnames = {'spprod','small phyto production','mg C m-3 day-1';
             'lpprod','large phyto production','mg C m-3 day-1';
             'bprod','bacteria production','mg C m-3 day-1';
             'szprod','small zooplankton production','mg C m-3 day-1';
             'mzprod','medium zooplankton production','mg C m-3 day-1';
             'lzprod','large zooplankton production','mg C m-3 day-1';
             'theta_sp','small phyto Chl:C','mg chl:mg C';
             'theta_lp','large phyto Chl:C','mg chl:mg C';
             'chl','chlorophyll','mg chl m-3'};
        
diag = zeros(size(Grd.z,1),size(diagnames,1));

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passsed variables are as follows:
% oldbio: array of state variable values to be updated
% par: photosynthetically available radiation for time step
% dmean_par: daily mean par for chl:C calculation
% kpar: attenuation coefficient for par
% temp = temperature
% z = depth coordinate
% dz = thickness of grid cell
% params = COBALT paramter values
% diag = diagnotics to be filled in
% t = time
% dt = time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = P.T;
sig = P.Sig;

nz = size(G.z,1);

sp    = oldbio(:,1);
lp    = oldbio(:,2);
sz    = oldbio(:,3);
mz    = oldbio(:,4);
lz    = oldbio(:,5);
b     = oldbio(:,6);
det   = oldbio(:,7);
ldon  = oldbio(:,8);
sldon = oldbio(:,9);
srdon = oldbio(:,10);
nh4   = oldbio(:,11);
no3   = oldbio(:,12);

newbio = zeros(size(oldbio));

% Rename all variables in the param structure to local variables
wsnk = B.wsnk;
ref_conc = B.ref_conc;
% Phytoplankton
kappa_phyto = B.kappa_phyto;
alpha_sp = B.alpha_sp;
alpha_lp = B.alpha_lp;
P_C_max_sp = B.P_C_max_sp;
P_C_max_lp = B.P_C_max_lp;
thetamax_sp = B.thetamax_sp;
thetamax_lp = B.thetamax_lp;
zeta_sp = B.zeta_sp;
zeta_lp = B.zeta_lp;
bresp_sp = B.bresp_sp;
bresp_lp = B.bresp_lp;
thetamin = B.thetamin;
kno3_sp = B.kno3_sp;
kno3_lp = B.kno3_lp;
knh4_sp = B.knh4_sp;
knh4_lp = B.knh4_lp;
m_agg_sp = B.m_agg_sp;
m_agg_lp = B.m_agg_lp;
m_vir_sp = B.m_vir_sp;
m_vir_lp = B.m_vir_lp;
exu = B.exu;
mld_thresh = B.mld_thresh;
% Zooplankton
kappa_zoo = B.kappa_zoo;
i_max_sz = B.i_max_sz;
i_max_mz = B.i_max_mz;
i_max_lz = B.i_max_lz;
ki_sz = B.ki_sz;
ki_mz = B.ki_mz;
ki_lz = B.ki_lz;
gge_max_sz = B.gge_max_sz;   
gge_max_mz = B.gge_max_mz;
gge_max_lz = B.gge_max_lz;          
bresp_sz = B.bresp_sz; 
bresp_mz = B.bresp_mz;
bresp_lz = B.bresp_lz;
phi_det_sz = B.phi_det_sz;
phi_det_mz = B.phi_det_mz;
phi_det_lz = B.phi_det_lz;
phi_ldon_sz = B.phi_ldon_sz;
phi_ldon_mz = B.phi_ldon_mz;
phi_ldon_lz = B.phi_ldon_lz;
phi_sldon_sz = B.phi_sldon_sz;
phi_sldon_mz = B.phi_sldon_mz;
phi_sldon_lz = B.phi_sldon_lz;
phi_srdon_sz = B.phi_srdon_sz;
phi_srdon_mz = B.phi_srdon_mz;
phi_srdon_lz = B.phi_srdon_lz;
ms_sz = B.ms_sz;
ns_sz = B.ns_sz;
ms_mz = B.ms_mz;
ns_mz = B.ns_mz;
ms_lz = B.ms_lz;
ns_lz = B.ns_lz;
ipa_matrix_sz = B.ipa_matrix_sz;
ipa_matrix_mz = B.ipa_matrix_mz;
ipa_matrix_lz = B.ipa_matrix_lz;
% Higher predation closure
kappa_hp = B.kappa_hp;
i_max_hp = B.i_max_hp;
ki_hp = B.ki_hp;
ms_hp = B.ms_hp;
ns_hp = B.ns_hp;
ipa_matrix_hp = B.ipa_matrix_hp;
phi_det_hp = B.phi_det_hp;
phi_ldon_hp = B.phi_ldon_hp;
phi_sldon_hp = B.phi_sldon_hp;
phi_srdon_hp = B.phi_srdon_hp;
phi_nh4_hp = B.phi_nh4_hp;
coef_hp = B.coef_hp;
% bacteria
kappa_bact = B.kappa_bact;
mu_max_b = B.mu_max_b;             
kldon = B.kldon;
gge_max_b = B.gge_max_b;
bresp_b = B.bresp_b;
m_vir_b = B.m_vir_b;
phi_ldon_vir = B.phi_ldon_vir;
phi_sldon_vir = B.phi_sldon_vir;
phi_srdon_vir = B.phi_srdon_vir; 
% Detritus and dissolved organic material
gamma_det = B.gamma_det;
gamma_srdon = B.gamma_srdon;
gamma_sldon = B.gamma_sldon;
gamma_nitrif = B.gamma_nitrif;

% Light within water column
parz = P.par.*exp(P.kpar.*G.z);    % depth-dependant irradiance (z is negative)
parz24 = P.par24.*exp(P.kpar.*G.z);  % irradiance integrated over 24 hours
% Chl:C based on 24 hour light exposure with light levels homogenized in
% the surface mixed layer and decaying exponentially below
sigdifcum = cumsum(diff(sig));      % cumulative dens diff from surface
aa = find(sigdifcum > mld_thresh);
parz24_mixed = parz24;
parz24_mixed(1:aa(1)) = mean(parz(1:aa(1)));
mld = aa(1)*G.dz;

% Calculate temperature scaling for each rate.  The parameter values for
% all rates are 0 deg. C.
tfac_sp = exp(kappa_phyto*temp);
tfac_lp = exp(kappa_phyto*temp);
tfac_sz = exp(kappa_zoo*temp);
tfac_mz = exp(kappa_zoo*temp);
tfac_lz = exp(kappa_zoo*temp);
tfac_hp = exp(kappa_hp*temp);
tfac_b =  exp(kappa_bact*temp);

% calculate refuge factors (applied to basal metabolic rates)
ref_sp = sp./(ref_conc + sp);
ref_lp = lp./(ref_conc + lp);
ref_sz = sz./(ref_conc + sz);
ref_mz = mz./(ref_conc + mz);
ref_lz = lz./(ref_conc + lz);
ref_b = b./(ref_conc + b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the phytoplankton growth rates                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the nutrient limitation (Frost and Franzen)
no3lim_sp = no3./( (kno3_sp + no3).*(1 + nh4/knh4_sp) );
no3lim_lp = no3./( (kno3_lp + no3).*(1 + nh4/knh4_lp) );
nh4lim_sp = nh4./(knh4_sp + nh4);
nh4lim_lp = nh4./(knh4_lp + nh4);
nlim_sp = no3lim_sp+nh4lim_sp;
nlim_lp = no3lim_lp+nh4lim_lp;

% calculate the maximum photosynthetic rate at a given temp/nutrients
P_C_m_sp = tfac_sp.*nlim_sp.*P_C_max_sp;
P_C_m_lp = tfac_lp.*nlim_lp.*P_C_max_lp;

% calculate the chlorophyll a:carbon ratio (modified Geider to include
% theta_min)
theta_sp = (thetamax_sp-thetamin)./ (1 + thetamax_sp.*alpha_sp.*parz24_mixed./ ...
                          (2.*P_C_m_sp)) + thetamin;
theta_lp = (thetamax_lp-thetamin)./ (1 + thetamax_lp.*alpha_lp.*parz24_mixed./ ...
                          (2.*P_C_m_lp)) + thetamin;

chl = theta_sp.*sp*6.625*12 + theta_lp.*lp*(106/16)*12;

% calculate the growth rates (Geider)
mu_sp = (P_C_m_sp./(1+zeta_sp) ).*(1-exp(-alpha_sp.*parz.*theta_sp./P_C_m_sp)) - ...
        tfac_sp.*bresp_sp.*ref_sp;
mu_lp = (P_C_m_lp./(1+zeta_lp) ).*(1-exp(-alpha_lp.*parz.*theta_lp./P_C_m_lp)) - ...
        tfac_lp.*bresp_lp.*ref_lp;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nutrient uptake calculations                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

juptake_no3_sp = max(0.0, mu_sp.*sp.*(no3lim_sp./(no3lim_sp + nh4lim_sp)));
juptake_no3_lp = max(0.0, mu_lp.*lp.*(no3lim_lp./(no3lim_lp + nh4lim_lp)));
juptake_nh4_sp = max(0.0, mu_sp.*sp.*(nh4lim_sp./(no3lim_sp + nh4lim_sp)));
juptake_nh4_lp = max(0.0, mu_lp.*lp.*(nh4lim_lp./(no3lim_lp + nh4lim_lp)));
% if growth is negative, create a negative juptake_nh4 equivalent to the
% total losses (i.e., net respiration and excretion occurs)
juptake_nh4_sp = juptake_nh4_sp + min(0.0,mu_sp.*sp);
juptake_nh4_lp = juptake_nh4_lp + min(0.0,mu_lp.*lp);

% production diagnostics
spprod = (juptake_no3_sp + juptake_nh4_sp);
lpprod = (juptake_no3_lp + juptake_nh4_lp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bacteria Dynamics                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate maximum bacterial uptake from specified maximum growth rate and
% maximum growth efficiency
vmax_b = (1.0/gge_max_b).*(mu_max_b + bresp_b);
bact_ldon_lim = ldon./(kldon + ldon);
juptake_ldon = vmax_b.*tfac_b.*bact_ldon_lim.*b;
bprod = gge_max_b.*juptake_ldon - ref_b.*tfac_b.*bresp_b.*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planktonic food web dynamics                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zooplankton Consumption

% 3 zooplankton types feeding on seven prey types [sp lp sz mz lz b det]
% "innate" refers to the innate availablity of the prey before adjusting
% for density-dependent switch as described in Stock et al. (2008),
% Bottom-up and top-down forcing in a simple size-structured plankton
% dynamics model, JMS, 74(1-2).  Dimension is nz x 7
innate_prey_sz = ipa_matrix_sz.*oldbio(:,1:7);
innate_prey_mz = ipa_matrix_mz.*oldbio(:,1:7);
innate_prey_lz = ipa_matrix_lz.*oldbio(:,1:7);
innate_prey_hp = ipa_matrix_hp.*oldbio(:,1:7);

% Create a vector with the denominators for the density dependent
% switching factors for each depth (nz x 1 column).  If ns = 1, this is
% just the sum of all available prey at each depth.
switch_norm_sz = sum(innate_prey_sz.^ns_sz,2);
switch_norm_mz = sum(innate_prey_mz.^ns_mz,2);
switch_norm_lz = sum(innate_prey_lz.^ns_lz,2);
switch_norm_hp = sum(innate_prey_hp.^ns_hp,2);

% make an nz x 7 matrix with the density dependent switching factors.  If
% ns = ms = 1 this reduces to the ratio of an individual prey type to the
% total prey and produces a very strong density-dependent switching
% response.  The default (ms = ns = 2) weakens this response considerably
switch_fac_sz = ipa_matrix_sz.* ...
                ( innate_prey_sz.^ns_sz./(switch_norm_sz*ones(1,7)) ).^(1/ms_sz);

switch_fac_mz = ipa_matrix_mz.* ...
                ( innate_prey_mz.^ns_mz./(switch_norm_mz*ones(1,7)) ).^(1/ms_mz);

switch_fac_lz = ipa_matrix_lz.* ...
                ( innate_prey_lz.^ns_lz./(switch_norm_lz*ones(1,7)) ).^(1/ms_lz);
            
switch_fac_hp = ipa_matrix_hp.* ...
                ( innate_prey_hp.^ns_hp./(switch_norm_hp*ones(1,7)) ).^(1/ms_hp);
            
% calculate the total prey available to each zooplankton group at each
% depth level (nz x 1 array)
available_prey_sz = sum(switch_fac_sz.*oldbio(:,1:7),2);
available_prey_mz = sum(switch_fac_mz.*oldbio(:,1:7),2);
available_prey_lz = sum(switch_fac_lz.*oldbio(:,1:7),2);
available_prey_hp = sum(switch_fac_hp.*oldbio(:,1:7),2);

% zooplankton ingestion rate matrix
ingrate_matrix_sz = (tfac_sz*ones(1,7)).*i_max_sz.*switch_fac_sz.*oldbio(:,1:7)./ ...
                    (ki_sz + available_prey_sz*ones(1,7));
ingrate_matrix_mz = (tfac_mz*ones(1,7)).*i_max_mz.*switch_fac_mz.*oldbio(:,1:7)./ ...
                    (ki_mz + available_prey_mz*ones(1,7));
ingrate_matrix_lz = (tfac_lz*ones(1,7)).*i_max_lz.*switch_fac_lz.*oldbio(:,1:7)./ ...
                    (ki_lz + available_prey_lz*ones(1,7));
ingrate_matrix_hp = (tfac_hp*ones(1,7)).*i_max_hp.*switch_fac_hp.*oldbio(:,1:7)./ ...
                    (ki_hp + available_prey_hp*ones(1,7));
                
% If there is no prey, calc above yield NaNs (Divide by 0 in switching
% factor).  Correct this by sumsetting the ingestion rate to 0 when no prey is
% present.
aa = find(isnan(ingrate_matrix_sz)); ingrate_matrix_sz(aa) = 0; clear aa;
aa = find(isnan(ingrate_matrix_mz)); ingrate_matrix_mz(aa) = 0; clear aa;
aa = find(isnan(ingrate_matrix_lz)); ingrate_matrix_lz(aa) = 0; clear aa;
aa = find(isnan(ingrate_matrix_hp)); ingrate_matrix_hp(aa) = 0; clear aa;

% total ingestion by each zooplankton type (nz x 1 column)
jingest_sz = sum(ingrate_matrix_sz,2).*sz;
jingest_mz = sum(ingrate_matrix_mz,2).*mz;
jingest_lz = sum(ingrate_matrix_lz,2).*lz;
% biomass of unresolved predator assumed proportional to the biomass of
% available prey (i.e., quadratic relationship, coef_hp = 2)
jingest_hp = sum(ingrate_matrix_hp,2).*available_prey_hp.*(coef_hp-1.0);

% losses to zooplankton ingestion
jzloss_sp = ingrate_matrix_sz(:,1).*sz + ingrate_matrix_mz(:,1).*mz + ...
            ingrate_matrix_lz(:,1).*lz;
jzloss_lp = ingrate_matrix_sz(:,2).*sz + ingrate_matrix_mz(:,2).*mz + ...
            ingrate_matrix_lz(:,2).*lz;
jzloss_sz = ingrate_matrix_sz(:,3).*sz + ingrate_matrix_mz(:,3).*mz + ...
            ingrate_matrix_lz(:,3).*lz;
jzloss_mz = ingrate_matrix_sz(:,4).*sz + ingrate_matrix_mz(:,4).*mz + ...
            ingrate_matrix_lz(:,4).*lz;
jzloss_lz = ingrate_matrix_sz(:,5).*sz + ingrate_matrix_mz(:,5).*mz + ...
            ingrate_matrix_lz(:,5).*lz;
jzloss_b = ingrate_matrix_sz(:,6).*sz + ingrate_matrix_mz(:,6).*mz + ...
           ingrate_matrix_lz(:,6).*lz;
jzloss_det = ingrate_matrix_sz(:,7).*sz + ingrate_matrix_mz(:,7).*mz + ...
             ingrate_matrix_lz(:,7).*lz;

% losses to higher predator ingestion
jhploss_sp = ingrate_matrix_hp(:,1).*available_prey_hp.*(coef_hp-1.0);
jhploss_lp = ingrate_matrix_hp(:,2).*available_prey_hp.*(coef_hp-1.0);
jhploss_sz = ingrate_matrix_hp(:,3).*available_prey_hp.*(coef_hp-1.0);
jhploss_mz = ingrate_matrix_hp(:,4).*available_prey_hp.*(coef_hp-1.0);
jhploss_lz = ingrate_matrix_hp(:,5).*available_prey_hp.*(coef_hp-1.0);
jhploss_b = ingrate_matrix_hp(:,6).*available_prey_hp.*(coef_hp-1.0);
jhploss_det = ingrate_matrix_hp(:,7).*available_prey_hp.*(coef_hp-1.0);

% losses to aggregation
jaggloss_sp = m_agg_sp.*sp.^2;
jaggloss_lp = m_agg_lp.*lp.^2;

% losses to viruses
jvirloss_sp = tfac_b*m_vir_sp.*sp.^2;
jvirloss_lp = tfac_b*m_vir_lp.*lp.^2;
jvirloss_b = tfac_b*m_vir_b.*b.^2;

% exudation
jexuloss_sp = exu.*max(juptake_no3_sp + juptake_nh4_sp,0.0);
jexuloss_lp = exu.*max(juptake_no3_lp + juptake_nh4_lp,0.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Production Calculations                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sinking detritus and dissolved organic material, note that negative
% bacterial production is allocated in the same proportions as viral 
% lysis.
jprod_det = phi_det_sz*jingest_sz + phi_det_mz*jingest_mz + ...
             phi_det_lz*jingest_lz + phi_det_hp*jingest_hp + ...
             jaggloss_sp + jaggloss_lp;
jprod_ldon = phi_ldon_sz*jingest_sz + phi_ldon_mz*jingest_mz + ...
             phi_ldon_lz*jingest_lz + phi_ldon_hp*jingest_hp + ...
             jexuloss_sp + jexuloss_lp + phi_ldon_vir*(jvirloss_sp + ...
             jvirloss_lp + jvirloss_b - min(bprod,0)) + gamma_sldon*sldon + ...
             gamma_srdon*srdon;
jprod_sldon = phi_sldon_sz*jingest_sz + phi_sldon_mz*jingest_mz + ...
              phi_sldon_lz*jingest_lz + phi_sldon_hp*jingest_hp + ...
              phi_sldon_vir*(jvirloss_sp + jvirloss_lp + jvirloss_b - ...
              min(bprod,0));
jprod_srdon = phi_srdon_sz*jingest_sz + phi_srdon_mz*jingest_mz + ...
              phi_srdon_lz*jingest_lz + phi_srdon_hp*jingest_hp + ...
              phi_srdon_vir*(jvirloss_sp + jvirloss_lp + jvirloss_b - ...
              min(bprod,0));

% start summing production terms for nh4 with remineralization of 
% dissolved organic material by bacteria
jprod_nh4 = juptake_ldon - max(bprod,0);

% Zooplankton production and excretion calculations, ingested material 
% that does not contribute to production or egestion is excreted, if 
% production is < 0, then negative production is allocated to sinking detritus.
% Small Zooplankton
szprod = gge_max_sz*jingest_sz - ref_sz.*tfac_sz.*bresp_sz.*sz;
aa = find(szprod > 0);
jprod_nh4(aa) = jprod_nh4(aa) + (jingest_sz(aa) - phi_det_sz*jingest_sz(aa) - ...
    phi_ldon_sz*jingest_sz(aa) - phi_sldon_sz*jingest_sz(aa) - ...
    phi_srdon_sz*jingest_sz(aa) - szprod(aa));
bb = find(szprod <= 0);
jprod_nh4(bb) = jprod_nh4(bb) + (jingest_sz(bb) - phi_det_sz*jingest_sz(bb) - ...
    phi_ldon_sz*jingest_sz(bb) - phi_sldon_sz*jingest_sz(bb) - ...
    phi_srdon_sz*jingest_sz(bb));
jprod_det(bb) = jprod_det(bb) - szprod(bb);
clear aa bb;
% Medium Zooplankton
mzprod = gge_max_mz*jingest_mz - ref_mz.*tfac_mz.*bresp_mz.*mz;
aa = find(mzprod > 0);
jprod_nh4(aa) = jprod_nh4(aa) + (jingest_mz(aa) - phi_det_mz*jingest_mz(aa) - ...
    phi_ldon_mz*jingest_mz(aa) - phi_sldon_mz*jingest_mz(aa) - ...
    phi_srdon_mz*jingest_mz(aa) - mzprod(aa));
bb = find(mzprod <= 0);
jprod_nh4(bb) = jprod_nh4(bb) + (jingest_mz(bb) - phi_det_mz*jingest_mz(bb) - ...
    phi_ldon_mz*jingest_mz(bb) - phi_sldon_mz*jingest_mz(bb) - ...
    phi_srdon_mz*jingest_mz(bb));
jprod_det(bb) = jprod_det(bb) - mzprod(bb);
clear aa bb;
% Large Zooplankton
lzprod = gge_max_lz*jingest_lz - ref_lz.*tfac_lz.*bresp_lz.*lz;
aa = find(lzprod > 0);
jprod_nh4(aa) = jprod_nh4(aa) + (jingest_lz(aa) - phi_det_lz*jingest_lz(aa) - ...
    phi_ldon_lz*jingest_lz(aa) - phi_sldon_lz*jingest_lz(aa) - ...
    phi_srdon_lz*jingest_lz(aa) - lzprod(aa));
bb = find(lzprod <= 0);
jprod_nh4(bb) = jprod_nh4(bb) + (jingest_lz(bb) - phi_det_lz*jingest_lz(bb) - ...
    phi_ldon_lz*jingest_lz(bb) - phi_sldon_lz*jingest_lz(bb) - ...
    phi_srdon_lz*jingest_lz(bb));
jprod_det(bb) = jprod_det(bb) - lzprod(bb);
clear aa bb;
% Add excretion by higher predators
jprod_nh4 = jprod_nh4 + phi_nh4_hp*jingest_hp;

% Remineralization of sinking detritus, mineral ballasting not yet
% implemented in matlab version
jremin_det = gamma_det*det;
jprod_nh4 = jprod_nh4 + jremin_det;

% nitrification (limitation of nitrification presumed to go as the small
% phytoplankton ammonia limitation)
jnitrif = gamma_nitrif*tfac_b.*nh4.*nh4lim_sp;

% Update state variables for biological sources and sinks
% Small Phytoplankton
newbio(:,1) = oldbio(:,1) + (spprod - jzloss_sp - jaggloss_sp - ...
    jvirloss_sp - jexuloss_sp - jhploss_sp).*G.dt;
% Large phytoplankton
newbio(:,2) = oldbio(:,2) + (lpprod - jzloss_lp - jaggloss_lp - ...
    jvirloss_lp - jexuloss_lp - jhploss_lp).*G.dt;
% Small Zooplankton
newbio(:,3)  = oldbio(:,3) + (szprod - jzloss_sz - jhploss_sz).*G.dt;
% Medium Zooplankton
newbio(:,4)  = oldbio(:,4) + (mzprod - jzloss_mz - jhploss_mz).*G.dt;
% Large Zooplankton
newbio(:,5)  = oldbio(:,5) + (lzprod - jzloss_lz - jhploss_lz).*G.dt;
% Bacteria
newbio(:,6) = oldbio(:,6) + (bprod - jzloss_b - jvirloss_b - jhploss_b).*G.dt;
% Sinking detritus
newbio(:,7) = oldbio(:,7) + (jprod_det - jzloss_det - jhploss_det - ...
    jremin_det).*G.dt;
% Labile dissolved organic material
newbio(:,8) = oldbio(:,8) + (jprod_ldon - juptake_ldon).*G.dt;
% Semi-labile dissolved organic material
newbio(:,9) = oldbio(:,9) + (jprod_sldon - gamma_sldon*sldon).*G.dt;
% Semi-refractory dissolved organic material
newbio(:,10) = oldbio(:,10) + (jprod_srdon - gamma_srdon*srdon).*G.dt;
% Ammonia
newbio(:,11) = oldbio(:,11) + (jprod_nh4 - juptake_nh4_sp - ...
               juptake_nh4_lp - jnitrif).*G.dt;
% Nitrate
newbio(:,12) = oldbio(:,12) + (jnitrif - juptake_no3_sp - juptake_no3_lp).*G.dt;

if any(newbio < 0)
    error('Negative tracer');
end

%sum(newbio(:))

diag = zeros(size(spprod,1),9);

diag(:,1) = spprod*86400*12*(106/16);
diag(:,2) = lpprod*86400*12*(106/16);
diag(:,3) = bprod*86400*12*(106/16);
diag(:,4) = szprod*86400*12*(106/16);
diag(:,5) = mzprod*86400*12*(106/16);
diag(:,6) = lzprod*86400*12*(106/16);
diag(:,7) = theta_sp;
diag(:,8) = theta_lp;
diag(:,9) = chl;

%diagnames = {'spprod','small phyto production','mg C m-3 day-1';
%             'lpprod','large phyto production','mg C m-3 day-1';
%             'bprod','bacteria production','mg C m-3 day-1';
%             'szprod','small zooplankton production','mg C m-3 day-1';
%             'mzprod','medium zooplankton production','mg C m-3 day-1';
%             'lzprod','large zooplankton production','mg C m-3 day-1';
%             'theta_sp','small phyto Chl:C','mg chl:mg C';
%             'theta_lp','large phyto Chl:C','mg chl:mg C';
%             'chl','chlorophyll','mg chl m-3'};


%**************************************************************************

function wsink = vertmove(bio, P, B, G)

% Creates a 2D field of vertical velocities for tracers.  Right now we just
% sink detritus, but more advanced parameterizations to handle, for
% example, vertical migration, are possible (e.g., Bianchi et al., 2013,
% GBC, 27  DOI:10.1002/gbc.20031)

wsink = zeros(size(bio));
wsink(:,7) = B.wsnk;