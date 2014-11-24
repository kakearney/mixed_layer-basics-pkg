function [db, Flx, Diag] = nemurokakode(time, bio, A)
%NEMUROKAKODE My version of NEMURO
%
% [db, Flx, Diag] = nemurokakode(time, bio, A)
%
% Source/sink ODE function for the nemurokak biological module.  See
% biomodules/nemurokak.m for details; this function is designed to be
% called by ODE solvers.


%------------------------------
% Set up various bio arrays
%------------------------------

[bv, ba, basum, bfrac, zlfrac, nb, nz] =  biomasssetup(bio, A);

%------------------------------
% Photosynthesis-related
% fluxes (gpp, exc, resp)
%------------------------------

[gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, felim, I, ...
    kappa, kappaP, fe2n, fedef] = primprod(bv.orig, A, nz, nb);

%------------------------------
% Grazing/predation
%------------------------------

graze = zeros(nb+2,nb+2,nz);
switch A.ivlev
    case 'single'
        graze(1:11,1:11,:) = ivlev(bv.orig(:,1:11), A.grmax, A.lambda, A.thresh, ...
                                    A.Kgra, A.temp);
    case 'multi'
        graze(1:11,1:11,:) = ivlevmulti(bv.orig(:,1:11), A.grmax, A.lambda, A.thresh, ...
                                    A.Kgra, A.temp, A.p);
    case 'orig'
        graze(1:11,1:11,:) = ivlevorignew(bv.prey(1:11,1:11,:), bv.pred(1:11,1:11,:), A.grmax, ...
                A.lambda, A.thresh, A.Kgra, A.temp, [A.psiPL A.psiZS]);
        if A.diapause

            % Calculate grazing using full population, and
            % hide-as-necessary population

%             [gratmp1, gratmp2] = deal(zeros(size(graze)));
            
            
%             gratmp1(1:11,1:11,:) = ivlevorig(bvin(:,1:11), A.grmax, ...
%                 A.lambda, A.thresh, A.Kgra, A.temp, [A.psiPL A.psiZS]);
%             gratmp2(1:11,1:11,:) = ivlevorig(bvgra(:,1:11), A.grmax, ...
%                 A.lambda, A.thresh, A.Kgra, A.temp, [A.psiPL A.psiZS]);
            
%             graze = gratmp1;
            
            % Then distribute both proportionally across the 2 groups
            
            zlfood = graze(:,A.idx.zl,:); 
            zlloss = graze(A.idx.zl,:,:); 
            
            zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [3 2 1]), [nb+2, 1, 1]));
            zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [2 3 1]), [1, nb+2, 1]));
            
            graze(:,[A.idx.zl1 A.idx.zl2],:) = zlfood;
            graze([A.idx.zl1 A.idx.zl2],:,:) = zlloss;    
            
            graze(:,A.idx.zl,:) = 0;
            graze(A.idx.zl,:,:) = 0;
            
%         else
%             graze(1:11,1:11,:) = ivlevorig(bv(:,1:11), A.grmax, A.lambda, A.thresh, ...
%                                         A.Kgra, A.temp, [A.psiPL A.psiZS]);
        end
    case 'mishmash'
        graze(1:11,1:11,:) = ivlevmishmash(bv(:,1:11), A.grmax, A.lambda, A.thresh, ...
                                    A.Kgra, A.temp, A.p, [A.psiPL A.psiZS]);
    
end

% No Si or Fe grazing fluxes.  Not tracking iron above phytoplankton, so
% grazing results in no net change.  Silica is all egested, so that will be
% dealt with below.  

pred = zeros(nb+2,nb+2,nz); % Just a placeholder for consistency w/ wce

%------------------------------
% Egestion and excretion by 
% consumers
%------------------------------

[egest, excrete, graze] =  egeexc(pred, graze, A, nb, nz, fe2n);

%------------------------------
% Mortality
%------------------------------

mort = nonpredmort(true, bv, A, fe2n, nb, nz);

%------------------------------
% Decomposition
%------------------------------

dec = decompremin(bv, A, nb, nz, I);

%------------------------------
% Gather all fluxes
%------------------------------

Flx.dec = dec;
Flx.ege = egest;
Flx.exx = exc;
Flx.exc = excrete;
Flx.gpp = gpp;
Flx.gra = graze;
Flx.mor = mort;
Flx.pre = pred;
Flx.res = resp;

% Reroute fluxes as indicated by user

if isfield(A, 'reroute') && ~isempty(A.reroute)
    
    for ir = 1:size(A.reroute,1)

        type = A.reroute{ir,1};
        idx1 = A.reroute{ir,2};
        idx2 = A.reroute{ir,3};
        idx3 = A.reroute{ir,4};
        frac = A.reroute{ir,5};

        reroute = Flx.(type)(idx1,idx2,:) .* frac;
        Flx.(type)(idx1,idx2,:) = Flx.(type)(idx1,idx2,:) - reroute;
        Flx.(type)(idx1,idx3,:) = Flx.(type)(idx1,idx3,:) + reroute;

    end
    
end

%------------------------------
% Total fluxes
%------------------------------  

fluxtot = Flx.dec + Flx.ege + Flx.exx + Flx.exc + Flx.gpp + ...
          Flx.gra + Flx.mor + Flx.pre + Flx.res;

fluxin  = permute(sum(fluxtot, 1), [3 2 1]);
fluxout = permute(sum(fluxtot, 2), [3 1 2]);

db = fluxin(:,1:end-2) - fluxout(:,1:end-2);

%------------------------------
% Diagnostics
%------------------------------ 

% If diapause flag is set, need to transfer the fluxes from the ZL1/ZL2
% groups to ZL

if A.diapause
    fx = fieldnames(Flx);
    for ii = 1:length(fx)
        Flx.(fx{ii})(A.idx.zl,:,:) = Flx.(fx{ii})(A.idx.zl1,:,:) + Flx.(fx{ii})(A.idx.zl2,:,:);
        Flx.(fx{ii})(:,A.idx.zl,:) = Flx.(fx{ii})(:,A.idx.zl1,:) + Flx.(fx{ii})(:,A.idx.zl2,:);
    end
end

% Phytoplankton growth limiters

Diag.psmax    = psmax( :, [A.idx.ps A.idx.pl]);
Diag.lightlim = Lfc(   :, [A.idx.ps A.idx.pl]);
Diag.no3lim   = no3lim(:, [A.idx.ps A.idx.pl]);
Diag.nh4lim   = nh4lim(:, [A.idx.ps A.idx.pl]);
Diag.silim    = silim( :,  A.idx.pl);
Diag.felim    = felim( :, [A.idx.ps A.idx.pl]);
Diag.templim  = bsxfun(@rdivide, psmax(:, [A.idx.ps A.idx.pl]), A.Vmax([A.idx.ps A.idx.pl])');

% Light-related stuff

Diag.I = I(:,1);
Diag.kappa = kappa;
Diag.kp = kappaP;

% Iron-related stuff

Diag.fe2n    = fe2n(   :, [A.idx.ps, A.idx.pl]);
Diag.fedef   = fedef(  :, [A.idx.ps A.idx.pl]);


%**************************************************************************

%------------------------
% Michaelis-Menten uptake
% limitation
%------------------------

function lim = mich(x, kx)
lim = x./(x + kx);
lim(isnan(lim)) = 0; % if x == 0 and kx == 0

function lim = mich2(x, kx)
lim = x.^2./(x.^2 + kx.^2);
lim(isnan(lim)) = 0; % if x == 0 and kx == 0

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);

%------------------------
% Single-resource Ivlev
%------------------------

function gra = ivlev(b, m, d, th, kgra, temp)
% d = nb x 1, predator 
% th = nb x nb, prey x pred
% b = nb x nz, both
% m = nb x nb, prey x pred
% kgra = nb x 1, predator
% temp = nz x 1

[nz,nb] = size(b);
[ipry,iprd] = find(m);

gra = zeros(nb,nb,nz);
for ii = 1:length(ipry)
    
    bprey = b(:,ipry(ii));
    bpred = b(:,iprd(ii));
    thresh = th(ipry(ii),iprd(ii));
    dd = d(iprd(ii));
    
    iv = max(1 - exp(dd .* (thresh - bprey)), 0);
    mm = tempdep(m(ipry(ii),iprd(ii)), kgra(iprd(ii)), temp);
    
    gra(ipry(ii),iprd(ii),:) = mm .* iv .* bpred;
end

%------------------------
% Multi-resource Ivlev
%------------------------

function gra = ivlevmulti(b, m, d, th, kgra, temp, p)
% d = nb x 1, predator 
% th = nb x 1, predator
% b = nb x nz, both
% m = nb x 1, pred
% p = nb x nb, prey x pred
% kgra = nb x 1, predator
% temp = nz x 1

[nz,nb] = size(b);
[ipry,iprd] = find(p);

gra = zeros(nb,nb,nz);
for ii = 1:length(ipry)
    
    weightedPrey = bsxfun(@times, b, p(:,iprd(ii))');
    sumPrey = sum(weightedPrey,2);
    fracPrey = bsxfun(@rdivide, weightedPrey, sumPrey);
    fracPrey = fracPrey(:,ipry(ii));
    
    bpred = b(:,iprd(ii));
    
    thresh = th(iprd(ii));
    dd = d(iprd(ii));
    
    iv = max((1 - exp(dd .* (thresh - sumPrey))).* fracPrey, 0);
    mm = tempdep(m(iprd(ii)), kgra(iprd(ii)), temp);
    
    gra(ipry(ii),iprd(ii),:) = mm .* iv .* bpred;
end

%------------------------
% Original Ivlev (single-
% resource for all but 
% ZP)
%------------------------

function gra = ivlevorig(b, m, d, th, kgra, temp, psi)
% d = nb x 1, predator 
% th = nb x nb, prey x pred
% b = nb x nz, both
% m = nb x nb, prey x pred
% kgra = nb x 1, predator
% temp = nz x 1
% psi = 2 x 1 (PL, ZS)

[nz,nb] = size(b);
[ipry,iprd] = find(m);

gra = zeros(nb,nb,nz);
for ii = 1:length(ipry)
    
    bprey = b(:,ipry(ii));
    bpred = b(:,iprd(ii));
    thresh = th(ipry(ii),iprd(ii));
    dd = d(iprd(ii));
    
    iv = max(1 - exp(dd .* (thresh - bprey)), 0);
    mm = tempdep(m(ipry(ii),iprd(ii)), kgra(iprd(ii)), temp);
    
    if iprd(ii) == 5 && ipry(ii) == 2 % PL -> ZP
        gourmet = exp(-psi(1) .* (b(:,3) + b(:,4)));
    elseif iprd(ii) == 5 && ipry(ii) == 3 % ZS -> ZP
        gourmet = exp(-psi(2) .* b(:,4));
    else
        gourmet = ones(nz,1);
    end
    
    gra(ipry(ii),iprd(ii),:) = mm .* iv .* gourmet .* bpred;
end

function gra = ivlevorignew(bi, bj, m, d, th, kgra, temp, psi)
% d = nb x 1, predator 
% th = nb x nb, prey x pred
% bi = nb x nb x nz, prey biomass
% bj = nb x nb x nz, predator biomass
% m = nb x nb, prey x pred
% kgra = nb x 1, predator
% temp = nz x 1
% psi = 2 x 1 (PL, ZS)

[nz,nb] = size(b);
[ipry,iprd] = find(m);

gra = zeros(nb,nb,nz);
for ii = 1:length(ipry)
    
    bprey = permute(bi(ipry(ii),iprd(ii),:), [3 1 2]);
    bpred = permute(bj(ipry(ii),iprd(ii),:), [3 1 2]);

%     bprey = b(:,ipry(ii));
%     bpred = b(:,iprd(ii));
    thresh = th(ipry(ii),iprd(ii));
    dd = d(iprd(ii));
    
    iv = max(1 - exp(dd .* (thresh - bprey)), 0);
    mm = tempdep(m(ipry(ii),iprd(ii)), kgra(iprd(ii)), temp);
    
    if iprd(ii) == 5 && ipry(ii) == 2 % PL -> ZP
        gourmet = exp(-psi(1) .* (b(:,3) + b(:,4)));
    elseif iprd(ii) == 5 && ipry(ii) == 3 % ZS -> ZP
        gourmet = exp(-psi(2) .* b(:,4));
    else
        gourmet = ones(nz,1);
    end
    
    gra(ipry(ii),iprd(ii),:) = mm .* iv .* gourmet .* bpred;
end

function gra = ivlevmishmash(b, m, d, th, kgra, temp, p, psi)
% d = nb x 1, predator 
% th = nb x 1, predator
% b = nb x nz, both
% m = nb x 1, pred
% p = nb x nb, prey x pred
% kgra = nb x 1, predator
% temp = nz x 1

[nz,nb] = size(b);
[ipry,iprd] = find(p);

gra = zeros(nb,nb,nz);
for ii = 1:length(ipry)
    
    weightedPrey = bsxfun(@times, b, p(:,iprd(ii))');
    sumPrey = sum(weightedPrey,2);
    fracPrey = bsxfun(@rdivide, weightedPrey, sumPrey);
    fracPrey = fracPrey(:,ipry(ii));
    
    bpred = b(:,iprd(ii));
    
    thresh = th(iprd(ii));
    dd = d(iprd(ii));
    
    iv = max((1 - exp(dd .* (thresh - sumPrey))).* fracPrey, 0);
    mm = tempdep(m(iprd(ii)), kgra(iprd(ii)), temp);
    
    % Use multi...
    
    gra(ipry(ii),iprd(ii),:) = mm .* iv .* bpred;
    
    % ...except for ZP
    
    bprey = b(:,ipry(ii));
    iv2 = max(1 - exp(dd .* (thresh - bprey)), 0);
    if iprd(ii) == 5 && ipry(ii) == 2 % PL -> ZP
        gourmet = exp(-psi(1) .* (b(:,3) + b(:,4)));
        gra(ipry(ii),iprd(ii),:) = mm .* iv2 .* bpred .* gourmet;
    elseif iprd(ii) == 5 && ipry(ii) == 3 % ZS -> ZP
        gourmet = exp(-psi(2) .* b(:,4));
        gra(ipry(ii),iprd(ii),:) = mm .* iv2 .* bpred .* gourmet;
    end
    
end
