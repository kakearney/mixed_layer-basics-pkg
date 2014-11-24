function [db, Flx, Diag] = wceode(time, bio, A)
%WCEODE Water column ecosystem model, main ODE function
%
% [db, Flx, Diag] = wceode(time, bio, A)
%
% Source/sink ODE function for wce module.  See biomodules/wce for details;
% this function is designed to be called by ODE solvers.


bio = max(bio, 0); % Negative biomass treated as 0

%------------------------------
% Set up various bio arrays
%------------------------------


% [bv,    bvin,    bvpre,    bvgra, ...
%  basum, basumin, basumpre, basumgra, ...
%  bfrac, bfracin, bfracpre, bfracgra, ...
%         zlfracin, zlfracpre, zlfracgra, nz, nb] = biomasssetup(bio, A);
    
[bv, ba, basum, bfrac, zlfrac, nb, nz] =  biomasssetup(bio, A);

%------------------------------
% Photosynthesis-related
% fluxes (gpp, exc, resp)
%------------------------------

[gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, felim, I, ...
    kappa, kappaP, fe2n, fedef] = primprod(bv.orig, A, nz, nb);


%------------------------------
% Aydin Ecosim primary 
% production
%------------------------------

if A.ecosimppflag  
    npp = zeros(nb+2,nb+2,nz);
    istoodeep = A.z < -100;
    for ib = [A.idx.ps A.idx.pl]
        h = exp(-2) + 1;  % TODO
        ybio = bv.orig(:,ib) ./ A.b0v(ib);
        npp(A.idx.mys,ib,:) = A.p0v(ib) .* (h.*ybio)./(h-1+ybio);
        npp(A.idx.mys,ib,istoodeep) = 0;
    end
end
   
%------------------------------
% Predation (involves nekton)
%------------------------------

pred = zeros(nb+2,nb+2,nz);

% Predation flux based on aydin ecosim functional response
% (ZL never predate, so don't have to worry about separating in/out here)

pred1 = aydinfrnew(basum.prey, basum.pred, A.b0, A.q0, A.x, A.d, A.theta); % mol N m^-2 s^-1
% pred1 = aydinfr(basumpre, A.b0, A.q0, A.x, A.d, A.theta); % mol N m^-2 s^-1

% For nekton-nekton links, this total flux sits a depth 1

prednn = zeros(nb);
prednn(A.links == 3) = pred1(A.links == 3);
prednn = prednn ./ dz(1); % mol N m^-3 s^-1

% For nekton-zooplankton links, the flux is distributed based on where the
% prey was located

prednz = bsxfun(@times, pred1, bfrac.prey);
prednz(isnan(prednz)) = 0;

% prednz = zeros(nb,nb,nz);
% for ipy = 1:nb
%     tmp1 = pred1(ipy,:);
%     tmp1(isnan(tmp1)) = 0;
%     tmp2 = bfrac.prey(:,ipy);
%     tmp2(isnan(tmp2)) = 0;
%     prednz(ipy,:,:) = (tmp2 * tmp1)'; 
% end

for iz = 1:nz
    prednz(:,:,iz) = prednz(:,:,iz)./dz(iz); % mol N m^-3 s^-1
end

for iprey = 1:nb
    for ipred = 1:nb
        if A.links(iprey,ipred) ~= 2
            prednz(iprey,ipred,:) = 0;
        end
    end
end

% Combine

pred(1:nb,1:nb,1) = prednz(:,:,1) + prednn;
pred(1:nb,1:nb,2:end) = prednz(:,:,2:end);

%------------------------------
% Grazing (involves only 
% plankton, and includes temp
% influence)
%------------------------------

% [graze1, graze2] = deal(zeros(nb+2,nb+2,nz));
% 
% for iz = 1:nz
%     graze1(1:nb,1:nb,iz) = aydinfrtemp(bvin(iz,:)',  A.b0v, A.q0vat0, A.x, A.d, A.theta, A.Kgra, A.temp(iz)); % mol N m^-3 s^-1
%     graze2(1:nb,1:nb,iz) = aydinfrtemp(bvgra(iz,:)', A.b0v, A.q0vat0, A.x, A.d, A.theta, A.Kgra, A.temp(iz)); % mol N m^-3 s^-1
% end

graze = zeros(nb+2,nb+2,nz);
for iz = 1:nz
    graze(1:nb,1:nb,iz) = aydinfrtempnew(bv.prey(:,:,iz), bv.pred(:,:,iz),  A.b0v, A.q0vat0, A.x, A.d, A.theta, A.Kgra, A.temp(iz)); % mol N m^-3 s^-1
end


for iprey = 1:nb
    for ipred = 1:nb
        if A.links(iprey,ipred) ~= 1 && A.links(iprey,ipred) ~= 4
%             graze1(iprey,ipred,:) = 0;
%             graze2(iprey,ipred,:) = 0;
            graze(iprey,ipred,:) = 0;
        end
    end
end

% graze = graze1;

% If diapause flag is on, split the ZL intake and loss proportionally
% across the ZL1 and ZL2 groups

if A.diapause
    
    % Grazing
    
%     zlfood = graze1(:,A.idx.zl,:);
%     zlloss = graze2(A.idx.zl,:,:);
    zlfood = graze(:,A.idx.zl,:);
    zlloss = graze(A.idx.zl,:,:);

%     zlfood = bsxfun(@times, zlfood, repmat(permute(zlfracin,  [3 2 1]), [nb+2, 1, 1]));
%     zlloss = bsxfun(@times, zlloss, repmat(permute(zlfracgra, [2 3 1]), [1, nb+2, 1]));
    zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [3 2 1]), [nb+2, 1, 1]));
    zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [2 3 1]), [1, nb+2, 1]));

    
    graze(:,[A.idx.zl1 A.idx.zl2],:) = zlfood;
    graze([A.idx.zl1 A.idx.zl2],:,:) = zlloss;    

    graze(:,A.idx.zl,:) = 0;
    graze(A.idx.zl,:,:) = 0;
    
    % Predation
    
    zlfood = pred(:,A.idx.zl,:);
    zlloss = pred(A.idx.zl,:,:);

%     zlfood = bsxfun(@times, zlfood, repmat(permute(zlfracpre, [3 2 1]), [nb+2, 1, 1]));
%     zlloss = bsxfun(@times, zlloss, repmat(permute(zlfracpre, [2 3 1]), [1, nb+2, 1]));
    zlfood = bsxfun(@times, zlfood, repmat(permute(zlfrac, [3 2 1]), [nb+2, 1, 1]));
    zlloss = bsxfun(@times, zlloss, repmat(permute(zlfrac, [2 3 1]), [1, nb+2, 1]));
    
    pred(:,[A.idx.zl1 A.idx.zl2],:) = zlfood;
    pred([A.idx.zl1 A.idx.zl2],:,:) = zlloss;    

    pred(:,A.idx.zl,:) = 0;
    pred(A.idx.zl,:,:) = 0;
    
end

%------------------------------
% Egestion and excretion
%------------------------------
    
[egest, excrete, graze] =  egeexc(pred, graze, A, nb, nz, fe2n);

%------------------------------
% Mortality
%------------------------------    

mort = nonpredmort(false, bv.orig, A, fe2n, nb, nz);

%------------------------------
% Decomposition
%------------------------------   

dec = decompremin(bv.orig, A, nb, nz, I);

%------------------------------
% Reroute
%------------------------------ 

% Gather together flux terms

if A.ecosimppflag
    Flx.dec = dec;
    Flx.ege = egest;
    Flx.exc = excrete;
    Flx.npp = npp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre = pred;
else
    Flx.dec = dec;
    Flx.ege = egest;
    Flx.exx = exc;
    Flx.exc = excrete;
    Flx.gpp = gpp;
    Flx.gra = graze;
    Flx.mor = mort;
    Flx.pre = pred;
    Flx.res = resp;;
end

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

if A.ecosimppflag
    fluxtot = Flx.dec + Flx.ege + Flx.exc + Flx.npp + ...
              Flx.gra + Flx.mor + Flx.pre;
else
    fluxtot = Flx.dec + Flx.ege + Flx.exx + Flx.exc + Flx.gpp + ...
              Flx.gra + Flx.mor + Flx.pre + Flx.res;
end
    
% fluxin  = squeeze(sum(fluxtot, 1))';
% fluxout = squeeze(sum(fluxtot, 2))';
fluxin  = permute(sum(fluxtot, 1), [3 2 1]);
fluxout = permute(sum(fluxtot, 2), [3 1 2]);

% Correct nekton fluxes so only enter and leave nekton in surface box
% (fluxes lower in the water column arise to grazing on and
% egesting/excreting plankton)

isnek = [A.isnek' false]; % add for critter Si

nekin = fluxin(:,isnek);
nekin = bsxfun(@times, nekin, dz);
nekin = sum(nekin,1)./dz(1);
fluxin(1,isnek) = nekin;
fluxin(2:end,isnek) = 0;

nekout = fluxout(:,isnek);
nekout = bsxfun(@times, nekout, dz);
nekout = sum(nekout,1)./dz(1);
fluxout(1,isnek) = nekout;
fluxout(2:end,isnek) = 0;

% Final rate of change

db = fluxin(:,1:nb) - fluxout(:,1:nb);
% db = fluxin - fluxout;

if any(isnan(db(:)))
    warning('WCE:nanInDbdt', 'NaN in dB/dt');
end


% if A.ecosimppflag
%     Flx.dec = dec;
%     Flx.egest = egest;
%     Flx.excrete = excrete;
%     Flx.npp = npp;
%     Flx.graze = graze;
%     Flx.mort = mort;
%     Flx.pred = pred;
% else
%     Flx.dec = dec;
%     Flx.egest = egest;
%     Flx.exc = exc;
%     Flx.excrete = excrete;
%     Flx.gpp = gpp;
%     Flx.graze = graze;
%     Flx.mort = mort;
%     Flx.pred = pred;
%     Flx.resp = resp;
% end

% If diapause flag is set, need to transfer the fluxes from the ZL1/ZL2
% groups to ZL

if A.diapause
    fx = fieldnames(Flx);
    for ii = 1:length(fx)
        Flx.(fx{ii})(A.idx.zl,:,:) = Flx.(fx{ii})(A.idx.zl1,:,:) + Flx.(fx{ii})(A.idx.zl2,:,:);
        Flx.(fx{ii})(:,A.idx.zl,:) = Flx.(fx{ii})(:,A.idx.zl1,:) + Flx.(fx{ii})(:,A.idx.zl2,:);
    end
end

% Diagnostics (all intermediate fluxes)

Diag.lightlim = Lfc(:,[A.idx.ps A.idx.pl]);
Diag.no3lim   = no3lim(:,[A.idx.ps A.idx.pl]);
Diag.nh4lim   = nh4lim(:,[A.idx.ps A.idx.pl]);
Diag.psmax    = psmax(:,[A.idx.ps A.idx.pl]);
Diag.silim    = silim(:,A.idx.pl);
Diag.I        = I(:,1);
Diag.kappa    = kappa;
Diag.kp       = kappaP;
Diag.felim    = felim(:,[A.idx.ps A.idx.pl]);
Diag.fe2n    = fe2n(  :,[A.idx.ps, A.idx.pl]);
Diag.fedef   = fedef( :,[A.idx.ps A.idx.pl]);

order = {'psmax', 'lightlim', 'no3lim', 'nh4lim', 'silim', 'felim', 'I', ...
         'kappa', 'kp', 'fe2n', 'fedef'};
Diag = orderfields(Diag, order);
     

flxcheck = structfun(@(x) any(x(:) < 0), Diag);
if any(flxcheck)
%     warning('WCE:negflux', 'Negative flux');
end

%------------------------------
% Aydin functional response
%------------------------------

function q = aydinfr(b, b0, q0, x, d, theta)
% Aydin-ecosim functional response
% b: nbsv x 1 vector
% units of b and q mut be consistent (either per area or per volume)!

nb = length(b);
b = b(:);

bi = b * ones(1,nb);
bj = ones(nb,1) * b';

b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0';

yi = bi./b0i;
yj = bj./b0j;

q = q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);

function q = aydinfrtemp(b, b0, q0, x, d, theta, kgra, temp)
% Aydin-ecosim functional response
% b: nbsv x 1 vector
% units of b and q mut be consistent (either per area or per volume)!

nb = length(b);
b = b(:);

bi = b * ones(1,nb);
bj = ones(nb,1) * b';

b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0';

kgra = ones(nb,1) * kgra';

yi = bi./b0i;
yj = bj./b0j;

q = exp(kgra.*temp) .* q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);

%------------------------------
% Aydin functional response
% with potentially different 
% prey and predator biomass
%------------------------------

function q = aydinfrnew(bi, bj, b0, q0, x, d, theta)

nb = size(bi,1);
b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0';

yi = bi./b0i;
yj = bj./b0j;

q = q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);


function q = aydinfrtempnew(bi, bj, b0, q0, x, d, theta, kgra, temp)

nb = size(bi,1);

b0i = b0 * ones(1,nb);
b0j = ones(nb,1) * b0';

kgra = ones(nb,1) * kgra';

yi = bi./b0i;
yj = bj./b0j;

q = exp(kgra.*temp) .* q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);


%------------------------------
% Aydin functional response
% (multi-dim version, was
% trying to speed things up but
% actually slower than the loop
% over each depth layer)
%------------------------------

function q = aydinfr2(b, b0, q0, x, d, theta)
% Aydin-ecosim functional response
% b:  nz x nbsv vector 
% b0: nbsv x 1
% q0: nbsv x nsbv vector
% units of b and q mut be consistent!

[nz,nb] = size(b);

yi = zeros(nb,nb,nz);
yj = zeros(nb,nb,nz);
for iz = 1:nz
    bi = b(iz,:)' * ones(1,nb);
    bj = ones(nb,1) * b(iz,:);

    b0i = b0 * ones(1,nb);
    b0j = ones(nb,1) * b0';
    
    yi(:,:,iz) = bi./b0i;
    yj(:,:,iz) = bj./b0j;
end

q0 = repmat(q0, [1 1 nz]);
x = repmat(x, [1 1 nz]);
d = repmat(d, [1 1 nz]);
theta = repmat(theta, [1 1 nz]);

q = q0 .* (x.*yj)./(x - 1 + yj) .* (d.*yi.^theta)./(d - 1 + yi.^theta);
    

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

    
    
    

