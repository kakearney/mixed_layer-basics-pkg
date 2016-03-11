function [gpp, exc, resp, dz, psmax, Lfc, no3lim, nh4lim, silim, felim, I, kappa, kappaP, fe2n, fedef] = primprod(bv, A, nz, nb)
%PRIMPROD Primary production fluxes for wce/nemurokak modules

if isscalar(A.dz)
    dz = ones(size(A.z)) * A.dz;
else
    dz = A.dz;
end

% Calculate light limitation factor

shadingbio = sum(bv(:,[A.idx.ps A.idx.pl]), 2);

nearsurf = 0.1*dz(1);    % Extend near to surface (since right at surface divides by 0)
zedge = [nearsurf; cumsum(dz)];
pintedge = [shadingbio(1).*nearsurf; cumsum(shadingbio).*dz];
kppedge = A.alpha2 .* pintedge./zedge;

kappaP = interp1(zedge, kppedge, -A.z);

kappa = A.alpha1 + kappaP;

I = A.irr .* exp(-kappa .* -A.z);

Lfc = zeros(nz,nb);
if A.usesteele
    I = I * ones(1,length(A.aIopt));
    Iopt = ones(nz,1) * A.Iopt';
    Lfc = I./Iopt .* exp(1 - I./Iopt);
else
    I = I * ones(1,length(A.alpha));
    alpha = ones(nz,1) * A.alpha';
    Vmax = ones(nz,1) * A.Vmax';
    Lfc = 1 - exp(-alpha .* I./Vmax);
end

% Calculate gross primary production

fratio = zeros(nz,nb);
[gpp, exc, resp] = deal(zeros(nb+2,nb+2,nz));

for ib = [A.idx.ps A.idx.pl]
    
    % Macronutrient limitation
    
    no3lim(:,ib) = mich(bv(:,A.idx.no3), A.Kno3(ib)) .* exp(-A.pusai(ib) .* bv(:,A.idx.nh4));
    nh4lim(:,ib) = mich(bv(:,A.idx.nh4), A.Knh4(ib));
    silim(:,ib)  = mich(bv(:,A.idx.sioh4), A.Ksi(ib)); % Always 1 if Ksi == 0, i.e. not si-dependant
%     felim(:,ib)  = mich(bv(:,A.idx.fe), A.Kfe(ib));
    
    % Iron limitation (quota model)
    
    if ib == A.idx.ps
        idxfe = A.idx.psfe;
    elseif ib == A.idx.pl
        idxfe = A.idx.plfe;
    end
    
%     fe2n(:,ib) = min(A.fe2nmax(ib), bv(:,idxfe)./bv(:,ib));  % molFe/molN
    fe2n(:,ib) = bv(:,idxfe)./bv(:,ib); % molFe/molN
    felim(:,ib) = mich(bv(:,A.idx.fe), A.Kfe(ib));           % Iron uptake limitation
    fedef(:,ib) = mich2(fe2n(:,ib), A.kfe2n(ib));            % Internal iron deficiency
    
    % F-ratio (fraction new production)
    
    fratio(:,ib) = no3lim(:,ib)./(no3lim(:,ib) + nh4lim(:,ib));
  
    % Overal nutrient limitation
    
%     nutlim = min([(no3lim(:,ib) + nh4lim(:,ib)) silim(:,ib) felim(:,ib)], [], 2);
    nutlim = min([(no3lim(:,ib) + nh4lim(:,ib)) silim(:,ib) fedef(:,ib)], [], 2);

    % Maximum possible production rate at current temps
    
    psmax(:,ib) = tempdep(A.Vmax(ib), A.Kgpp(ib), A.temp);
    
    % Nitrogen uptake
    
    totnuptake = psmax(:,ib) .* Lfc(:,ib) .* nutlim .* bv(:,ib);
    
    no3uptake = totnuptake .* fratio(:,ib);
    nh4uptake = totnuptake .* (1 - fratio(:,ib));
    
    % Uptake of SiOH4
    
    if ib == A.idx.pl
        siuptake = (no3uptake + nh4uptake) .* A.RSiN;
    else
        siuptake = zeros(nz,1);
    end
    
    % Uptake of Fe (quota model)
    
    isless = fe2n(:,ib) < A.fe2nmax(ib);
    feuptake(isless) = psmax(isless,ib) .* bv(isless,ib) .* felim(isless,ib) .* A.fe2nupfac;
    feuptake(~isless) = 0;
    
%     feuptake = (no3uptake + nh4uptake) .* (A.RCN/A.RCFe);
    
    % Gross primary production fluxes
    
    gpp(A.idx.no3,ib,:) = no3uptake;
    gpp(A.idx.nh4,ib,:) = nh4uptake;
    if ib == A.idx.pl
        gpp(A.idx.sioh4,A.idx.plsi,:) = siuptake;
    end
%     gpp(A.idx.fe,nb+1,:) = permute(feuptake, [2 3 1]) + gpp(A.idx.fe,folder+1,:);
    gpp(A.idx.fe,idxfe,:) = feuptake;

    % Extracellular excretion (N) is proportional to uptake
    
    excn = A.gamma(ib) .* (no3uptake + nh4uptake);
    exc(ib,A.idx.don,:) = excn;
        
    % Phytoplankton respiration (N)
    
    respn = tempdep(A.res0(ib), A.Kres(ib), A.temp) .* bv(:,ib);
    
    resp(ib,A.idx.no3,:) = respn .* fratio(:,ib);
    resp(ib,A.idx.nh4,:) = respn .* (1 - fratio(:,ib));
    
    % Excretion and respiration of silica proportional to N
    
    if ib == A.idx.pl
        excsi = excn * A.RSiN;
        exc(A.idx.plsi,A.idx.sioh4,:) = excsi;
        
        respsi = respn * A.RSiN;
        resp(A.idx.plsi, A.idx.sioh4,:) = respsi;
    end
    
    % Excretion and respiration of iron proportional to N
    % TODO: with quota model, should Fe be excreted and respired?
    
%     excfe = excn .* (A.RCN/A.RCFe);
%     exc(nb+1,A.idx.fe,:) = permute(excfe, [2 3 1]) + exc(nb+1,A.idx.fe,:);
%     
%     respfe = respn .* (A.RCN/A.RCFe);
%     resp(nb+1,A.idx.fe,:) = permute(respfe, [2 3 1]) + resp(nb+1,A.idx.fe,:);
    
end

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
