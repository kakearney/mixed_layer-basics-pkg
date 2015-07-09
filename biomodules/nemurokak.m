function varargout = nemurokak(action, varargin)
%NEMUROKAK NEMURO, with iron and grazing modifications
%
% See biomodule.m for full syntax details.
%
% This module runs a lower-trophic level biogeochemical model based on the
% NEMURO model.  It includes modifications for explicit iron limitation, as
% well as options to switch between single-resource and multi-resource
% grazing functional responses.
%
% See biomodule.m for function syntax descriptions.  For fields that must be 
% present in the In structure (passed to mixed_layer as parameter/value pairs), 
% see the help for parsewcenemin.m. 

% Copyright 2011-2015 Kelly Kearney

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

function [bio, ismixed, bottomval, Biovars, names, diagnames] = init(In, Grd)

%------------------------------
% Parse and check input
%------------------------------

BioIn = parsewcenemin('isnem', true, In);

%----------------------
% Variable names
%----------------------

[names, nbsv, nemidx, Biovars] = setstatevars(BioIn); 

[diagnames, Biovars] = setdiagnosticsvars(BioIn, names, Biovars);

%------------------------------
% Initial biomass
%------------------------------

S = warning('off', 'MATLAB:interp1:NaNinY');

% Check that bnem0 has correct # of columns

% nbsv = 18; % Number of tracked biological state variables

nz = length(Grd.z);
nb0 = size(BioIn.bnem0,2);
if nb0 ~= 13
    error('Initial biomass (bnem0) should include depth plus 12 state variables');
end

bio = zeros(nz, nbsv);

% Most state variables' initial profiles set by input...

if all(BioIn.bnem0 >= 0)
    BioIn.bnem0(:,1) = -BioIn.bnem0(:,1); % Tired of trying to remember +up or down
end

bio(:,1:12) = interp1(BioIn.bnem0(:,1), BioIn.bnem0(:,2:end), Grd.z); % mol N/m^3

% ... except PL silica, which is proportional to N

bio(:,13) = bio(:,2) .* BioIn.NemParam.RSiN;

% ... and PS and PL iron, which starts at 0.  TODO: Or should I assume
% something else?  Half-sat Fe:N maybe?

% bio(:,14:15) = bsxfun(@times, bio(:,1:2), BioIn.fe2nmax);
bio(:,14:15) = 0;

% ... and POFe starts at 0,

bio(:,16) = 0;

% If diapause option is on, split ZL into two groups.  Start with all
% biomass in the non-migrating group.

if BioIn.diapause
    zltot = bio(:,4);
    bio(:,17) = zltot;
    bio(:,18) = 0;
    bio(:,4)  = 0;
end

warning(S);

%----------------------
% Mixing and forcing
%----------------------

% All variables mixed

ismixed = true(1, nbsv);

% No bottom forcing

bottomval = nan(1, nbsv);

%------------------------------
% Variables needed for ODE
%------------------------------

% Params shared with wce

[Biovars, Np] = setwcenemparams(BioIn, nemidx, nbsv, Grd, Biovars);

% Specific to nemurokak

Biovars.grmax                             = zeros(nbsv,nbsv);
Biovars.grmax(nemidx(1:11),nemidx(1:11))  = Np.grmax;               % s^-1
Biovars.lambda                            = zeros(nbsv,1);
Biovars.lambda(nemidx(1:11))              = Np.lambda/1000;         % m^3/molN
Biovars.thresh                            = zeros(nbsv,nbsv);
Biovars.thresh(nemidx(1:11),nemidx(1:11)) = Np.thresh * 1000;       % molN/m^3
Biovars.Kgra                              = zeros(nbsv,1);
Biovars.Kgra(nemidx(1:11))                = Np.Kgra;                % degC^-1
Biovars.mor0                              = zeros(nbsv,1);
Biovars.mor0(nemidx(1:11))                = Np.mor0/1000;           % (molN/m^3)^-1 s^-1
Biovars.Kmor                              = zeros(nbsv,1);
Biovars.Kmor(nemidx(1:11))                = Np.Kmor;                % degC^-1

switch BioIn.ivlev

    case 'orig'

        % Need these ones for originl Ivlev only

        Biovars.psiPL = BioIn.NemParam.PusaiPL/1000; % m^3/molN
        Biovars.psiZS = BioIn.NemParam.PusaiZS/1000; % m^3/molN

    case 'multi' % Not up to date since diapause rewrite

        % User-input matrices overshadow NP set

        if isempty(BioIn.grmax) || isempty(BioIn.thresh) || isempty(BioIn.p)
            error('Need to provide grmax, thresh, and p for multi-resource');
        end

        Biovars.grmax  = BioIn.grmax;
        Biovars.thresh = BioIn.thresh;
        Biovars.p      = BioIn.p;

    case 'mishmash' % Not up to date since diapause rewrite

        Biovars.psiPL = BioIn.NemParam.PusaiPL/1000; % m^3/molN
        Biovars.psiZS = BioIn.NemParam.PusaiZS/1000; % m^3/molN
        if isempty(BioIn.grmax) || isempty(BioIn.thresh) || isempty(BioIn.p)
            error('Need to provide grmax, thresh, and p for multi-resource');
        end

        Biovars.grmax  = BioIn.grmax;
        Biovars.thresh = BioIn.thresh;
        Biovars.p      = BioIn.p;

end

Biovars.ivlev            = BioIn.ivlev;   
Biovars.gs               = zeros(nbsv,1);
Biovars.gs(nemidx(1:11)) = 1 - Np.alphaeg;   % no unit
Biovars.ge               = zeros(nbsv,1);
Biovars.ge(nemidx(1:11)) = Np.beta;          % no unit

if Biovars.diapause
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
                                      
% For diapause splitting/combining, need to do this outside the ODE

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

% Run main ODE

[newbio, db, Flx, Diag, badthings] = integratebio(@nemurokakode, G.t, G.dt, oldbio, Param, B.odesolver{:});

% Recombine ZL groups for physical mixing

if B.diapause
    newbio(:,B.idx.zl) = newbio(:,B.idx.zl1) + newbio(:,B.idx.zl2);
end

% Check for errors

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
    error('NEMURO:biologyOutOfRange', errstr);

end

% Diagnostics: Other

ndiag = 14;

if isempty(Diag) % If use any solver other than euler
    diag = zeros(length(z), ndiag);
else
    diag = struct2cell(Diag);
    diag = cat(2, diag{:});
end

% Diagnostics: Intermediate fluxes

nemflux = B.flux;
nfx = size(nemflux,1);
nz = size(newbio,1);
fluxes = zeros(nz, nfx);
for ifx = 1:nfx
    fluxes(:,ifx) = Flx(1).(nemflux{ifx,1})(nemflux{ifx,2}, nemflux{ifx,3},:);
end

diag = [diag fluxes];

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




