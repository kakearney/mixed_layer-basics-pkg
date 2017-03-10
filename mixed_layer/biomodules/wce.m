function varargout = wce(action, varargin)
%WCE Water column ecosystem biological module
%
% This module simulates a mixed planktonic-nektonic ecosystem, based on
% a combination of models derived mainly from NEMURO and Kerim Aydin's
% version of Ecosim, with a little bit of COBALT thrown in for flavor.
%
% A note on units: The biomass of all critters is saved to file in mol
% N[Si][Fe]/m^3.  For nektonic critters, all biomass is placed in the
% surface cell, and actually represents the total over the entire water
% column; multiply by the thickness of the surface layer to get the true
% biomass, in mol N/m^2.
%
% See biomodule.m for function syntax descriptions.  For fields that must be 
% present in the In structure (passed to mixed_layer as parameter/value pairs), 
% see the help for parsewcenemin.m. 

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

%------------------------------
% Parse and check input
%------------------------------

BioIn = parsewcenemin('isnem', false, In);

%------------------------------
% Variable names
%------------------------------

[names, nbsv, nemidx, Biovars] = setstatevars(BioIn); 

[diagnames, Biovars] = setdiagnosticsvars(BioIn, names, Biovars);

%----------------------------
% Initial biomass for
% state variables
%----------------------------

nz = length(Grd.z);

bio = zeros(nz, nbsv);

% Nemuro-derived non-live variables start with user-defined values in
% mol/m^3 (includes NO3 through Fe)

if all(BioIn.bnem0 >= 0)
    BioIn.bnem0(:,1) = -BioIn.bnem0(:,1); % Tired of trying to remember +up or down
end

bio(:,nemidx(6:12)) = interp1(BioIn.bnem0(:,1), BioIn.bnem0(:,7:13), Grd.z); % mol/m^3

% Plankton groups distribute their biomass evenly throughout the mixed
% layer (mld as defined by the user, not T and S)

isabovemld = Grd.z >= BioIn.mld;
nlayer = sum(isabovemld);
mld = -Grd.zp(nlayer+1);

EM = BioIn.EM;
if isempty(BioIn.ensdata)
    Ep = EM.ecopath;
else
    [~, Ep] = EM.ecopath('ensemble', BioIn.ensdata);
end

% Ewein = ecopathinputcheck(BioIn.Ewein, true);
% Ep = ecopathlite(Ewein);

isplank = Biovars.iszoo | Biovars.isphy;

bplank = Ep.b(isplank)'./mld;   % mol N/m^3
bio(isabovemld,isplank) = repmat(bplank, nlayer, 1); 

% Nekton biomass is stored in top layer for convenience.  Although actually
% per area, here stored as per volume for consistency

bnek = Ep.b(Biovars.isnek)'./-Grd.zp(2); % mol N/m^3
bio(1,Biovars.isnek) = bnek;

% PL silica is proportional to PL N

pln = bio(:, Biovars.idx.pl);
plsi = pln .* BioIn.NemParam.RSiN;
bio(:, Biovars.idx.plsi) = plsi;

% PS and PL iron starts at 0

bio(:,[Biovars.idx.psfe Biovars.idx.plfe]) = 0;

% So does POFe

bio(:, Biovars.idx.pofe) = 0;

% If diapause, split ZL biomass

if Biovars.diapause
    bio(:, Biovars.idx.zl1) = bio(:, Biovars.idx.zl);
    bio(:, Biovars.idx.zl) = 0;
end
    
%----------------------------
% Indicators for mixing and 
% bottom-forcing
%----------------------------

% Plankton mixed, nekton not
% If no-mix (debugging), nothing is mixed.

if BioIn.nomix
    ismixed = false(nbsv,1);
else
    ismixed = ~Biovars.isnek;
end

if Biovars.diapause
    ismixed([Biovars.idx.zl1 Biovars.idx.zl2]) = true;
end

% No bottom forcing

bottomval = nan(nbsv,1);

%----------------------------
% Variables needed for ODE
%----------------------------

% Params shared with nemurokak

[Biovars, Np] = setwcenemparams(BioIn, nemidx, nbsv, Grd, Biovars);

% Mortality exponent

Biovars.m0exp                = zeros(nbsv,1);
Biovars.m0exp(1:EM.nlive)    = BioIn.m0exp;            % no unit

% Nekton: the paramters for nekton come from the Ecopath mass balance

Biovars.b0 = zeros(nbsv,1);
Biovars.b0(1:EM.ngroup) = Ep.b;
Biovars.q0 = zeros(nbsv);
Biovars.q0(1:EM.nlive,1:EM.nlive) = Ep.q0(1:EM.nlive,1:EM.nlive);
Biovars.ge = zeros(nbsv,1);
Biovars.ge(1:EM.nlive) = Ep.ge(1:EM.nlive);
Biovars.gs = zeros(nbsv,1);
Biovars.gs(1:EM.nlive) = EM.groupdata.gs(1:EM.nlive);

% For x, d, and theta, I accept data in one of two ways:
% 1) vector of P-value log-tranformed-anomaly-from-base values, identical
% to those used in the functional response file input for aydin-ecosim
% 2) matrix of values for each group pair.  These values are NOT anomalies
% but actual values to be used in the functional response equations.

dcmask = table2array(EM.dc) == 0;

if isvector(BioIn.x)
    [xj,xi] = meshgrid(BioIn.x);
    x = 1 + exp(xi + xj);
    x(dcmask) = 0;
else
    x = BioIn.x;
end

if isvector(BioIn.d)
    [dj,di] = meshgrid(BioIn.d);
    d = 1 + exp(di + dj);
    d(dcmask) = 0;
else
    d = BioIn.d;
end

if isvector(BioIn.theta)
    [thj,thi] = meshgrid(BioIn.theta);
    theta = exp(0.05*(thi + thj)); % TODO double-check theta, since different in spreadsheet and ppt
    theta(dcmask) = 0;
else
    theta = BioIn.theta;
end

n = nbsv - EM.ngroup;
Biovars.x = padarray(x, [n n], 'post');
Biovars.d = padarray(d, [n n], 'post');
Biovars.theta = padarray(theta, [n n], 'post');

% Zooplankton: grazing same as predation, but in volumetric terms

dz = mld; % TODO: This is arbitrary, and possibly very important... run some tests

Biovars.b0v = Biovars.b0./dz; % mol N m^-3
Biovars.q0v = Biovars.q0./dz; % mol N m^-3 s^-1

% Production (for debugging with ecosim production only)

Biovars.p0 = zeros(nbsv,1);
Biovars.p0(1:EM.ngroup) = Ep.pb .* Ep.b; % mol N m^-2 s^-1
Biovars.p0v = Biovars.p0./dz;                     % mol N m^-3 s^-1 

% Temperature factors

Biovars.Kgra = zeros(nbsv,1);
Biovars.Kgra(Biovars.iszoo & ~Biovars.isextrazoo) = Np.Kgra(3:5);
Biovars.Kgra(Biovars.isextrazoo) = BioIn.kgra;

tempfac0 = exp(Biovars.Kgra .* BioIn.temp);

Biovars.q0vat0 = bsxfun(@rdivide, Biovars.q0v, tempfac0');

% Mortality (using ecopath-based mortality)

Biovars.m0 = zeros(nbsv,1);
Biovars.m0(1:EM.nlive) = Ep.otherMortRate(1:EM.nlive);  % s^-1
if any(Biovars.m0 < 0)
    warning('WCE:m0neg', 'Some m0 < 0, setting to 0');
    Biovars.m0 = max(Biovars.m0, 0);
end

% Mortality: Allows for any exponential function aB.^y

m0 = zeros(nbsv,1);
m0(1:EM.nlive) = Ep.otherMortRate(1:EM.nlive);  % s^-1
m0b = m0 .* bio(1,:)'; % mass-balanced flux, molN/m^3/s

Biovars.m0coef = m0b./(bio(1,:)'.^Biovars.m0exp);
Biovars.m0coef = max(Biovars.m0coef, 0); % get rid of NaNs (and I guess negatives, if they ever appeared, though they shouldn't)

if Biovars.diapause
    bzl = sum(bio(1,[Biovars.idx.zl1 Biovars.idx.zl2]));
    m0bzl = Ep.otherMortRate(Biovars.idx.zl) .* bzl;
    mexp = Biovars.m0exp(Biovars.idx.zl);
    Biovars.m0coef(Biovars.idx.zl) = 0;
    Biovars.m0coef([Biovars.idx.zl1 Biovars.idx.zl2]) = m0bzl./(bzl.^mexp);
    
    Biovars = setdiapauseparams(BioIn, Biovars, Grd);
end



%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)

% Set up parameters that will be passed to main ODE function

Param      = B;
Param.irr  = P.par24;
Param.temp = P.T;
Param.z    = G.z;
Param.dz   = G.dz;
    
if any(isnan(oldbio(:)))
    warning('WCE: NaN in biology');
end

% Split/combine ZL groups as necessary for diapause

if B.diapause
    it = find(G.t == B.t);
    zltot = sum(oldbio(:,[B.idx.zl1 B.idx.zl2]), 2);
    if B.zlsplit(it)
        ztransfer = oldbio(:, B.idx.zl1) * B.zlsplit(it);
        oldbio(:,B.idx.zl2) = oldbio(:,B.idx.zl2) + ztransfer;
        oldbio(:,B.idx.zl1) = oldbio(:,B.idx.zl1) - ztransfer;
    elseif B.zlcombine(it)
        oldbio(:,B.idx.zl1) = zltot;
        oldbio(:,B.idx.zl2) = 0;
    end
    oldbio(:,B.idx.zl) = 0;
end

% Integrate biology over this time step

[newbio, db, Flx, Diag, badthings] = integratebio(@wceode, G.t, G.dt, oldbio, Param, B.odesolver{:});

if B.diapause
    newbio(:,B.idx.zl) = newbio(:,B.idx.zl1) + newbio(:,B.idx.zl2);
end

% Check and correct for silica issue

isneg = newbio(:, B.idx.plsi) < 0;
Diag.extrasi = zeros(size(isneg));
Diag.extrasi(isneg) = -newbio(isneg,B.idx.plsi);
newbio(isneg, B.idx.plsi) = 0;
badthings(isneg, B.idx.plsi) = false;

% Check for problems

if any(badthings(:))
    [ridx,cidx] = find(badthings);
    nb = length(ridx);
    errstr = cell(nb,1);
    for ii = 1:nb
        badbio = newbio(ridx(ii), cidx(ii));
        if isnan(badbio)
            errstr{ii} = sprintf('NaN: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), G.t);
        elseif isinf(badbio)
            errstr{ii} = sprintf('Inf: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), G.t);
        elseif badbio < 0
            errstr{ii} = sprintf('Neg: depth %d, critter %d, time = %d', ridx(ii), cidx(ii), G.t); 
        end
    end
    errstr = sprintf('  %s\n', errstr{:});
    errstr = sprintf('Biology out of range:\n%s', errstr);
    error('WCE:biologyOutOfRange', errstr);
end

% Diagnostics: Intermediate fluxes

nd = size(B.flux,1);

fluxes = zeros(length(G.z),nd);
for id = 1:nd
    fluxes(:,id) = Flx(1).(B.flux{id,1})(B.flux{id,2},B.flux{id,3},:);
end

% Diagnistics: Other

ndiag = 15;
if isempty(Diag) % If use any solver other than euler
    diag = zeros(length(G.z), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(2, diag{:});
end

diag = [db fluxes diag];

% Conservation check (debugging, very ESA-specific)

% tot = sum(newbio .* Param.dz, 1); 
% ntot = sum(tot(1:27));          % mol m^-2
% stot = sum(tot([28:29 31]));    % mol m^-2
% ftot = sum(tot([30 32:33]));    % umol m^-2
% fprintf(Biovars.cfid, '%f %f %f %f\n', t, ntot, stot, ftot);
% if abs(ntot - Biovars.ntot0) > 0.01
%     fclose(Biovars.cfid);
%     error('WCE:notConserved', 'N not conserved');
% end
                      
%**************************************************************************

function wsink = vertmove(oldbio, P, B, G)

wsink = B.settle;

if B.diapause
    
    it = find(B.t == G.t);
    
    if B.zlswim(it) == 1
        wsink(G.z < -10, B.idx.zl2) = 80./86400; % swim up
    elseif B.zlswim(it) == -1
        wsink(G.z > -400, B.idx.zl2) = -80./86400; % stay down
        wsink(G.z < -450, B.idx.zl2) = 80./86400;
    end
    
end

