function dec = decompremin(bv, A, nb, nz, I)
%DECOMPREMIN Decomposition, remineralization, etc., for namurokak/wce

% NEMURO-derived decomposition/remineralization

dec = zeros(nb+2,nb+2,nz);
for isrc = 1:nb
    for isnk = 1:nb
        if A.vdec(isrc,isnk) > 0
            td = tempdep(A.vdec(isrc,isnk), A.Kdec(isrc,isnk), A.temp);
            dec(isrc,isnk,:) = td .* bv(:,isrc);          
        end
    end
end

% Fe scavenging

klig = min(A.kliglo, 10.^(log10(A.klighi) + max(0, log10(A.iofescav./I(:,1)))));

qa = klig;
qb = 1 + klig .* (A.ligbkg - bv(:,A.idx.fe));
qc = -bv(:,A.idx.fe);
fefree = (-qb + sqrt(qb.^2 - 4.*qa.*qc))./(2.*qa);  % Amount not bound to ligands

feads = (A.scavfac.*max(fefree-A.scavthresh,0) + A.alphascav) .* fefree;

% feads = A.alphascav .* fefree; % Adsorption onto particles of free iron

% TODO: What about when iron is high?  Higher scavenging?  When does
% precipitation start to occur?

dec(A.idx.fe, A.idx.pofe,:) = feads;


% POFe remineralization

td = tempdep(A.vdec(A.idx.pon,A.idx.nh4), A.Kdec(A.idx.pon,A.idx.nh4), A.temp);
dec(A.idx.pofe,A.idx.fe,:) = td .* bv(:,A.idx.pofe) .* A.remineff;

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);
