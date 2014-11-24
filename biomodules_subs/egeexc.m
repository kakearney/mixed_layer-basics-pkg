function [egest, excrete, graze] = egeexc(pred, graze, A, nb, nz, fe2n)

alleat = pred + graze; % (never overlap) nb+2 x nb+2 x nz
alleat = permute(sum(alleat,1), [3 2 1]);  %  nz x nb, total eaten by each predator

[egest, excrete] = deal(zeros(nb+2,nb+2,nz));

% Egestion and excretion

egetmp = bsxfun(@times, alleat, [A.gs' 0 0]);
exctmp = bsxfun(@times, alleat, [(1 - A.gs' - A.ge') 0 0]);

egest(:,A.idx.pon,:) = egetmp';
excrete(:,A.idx.nh4,:) = exctmp';
    
% Assume no critters can assimilate Si, all egested immediately to opal

grasi = graze(A.idx.pl,:,:) + pred(A.idx.pl,:,:); % 1 x nb+2 x nz
grasi = permute(grasi, [3 2 1]);                  %  nz x nb+2 
grasi = sum(grasi,2) .* A.RSiN;

egest(A.idx.plsi,A.idx.opal,:) = grasi;

% Assume iron assimilated, excreted, and egested proportional to N

% egefe = sum(egetmp, 2) .* (A.RCN/A.RCFe);
% egest(nb+1,A.idx.fe,:) = egefe;
% 
% excretefe = sum(exctmp, 2) .* (A.RCN/A.RCFe);
% excrete(nb+1,A.idx.fe,:) = excretefe;

% Assume no critter assimilate Fe, all egested to Fe

% grafe = graze([A.idx.ps A.idx.pl],:,:) + pred([A.idx.ps A.idx.pl],:,:);
% grafe = bsxfun(@times, grafe, permute(fe2n, [2 3 1]));
% grafe = sum(grafe, 2); 
% egest([A.idx.psfe A.idx.plfe], A.idx.fe, :) = grafe;

% Track egestion only as it passes to, and is egested/excreted by,
% first-level consumers.  Beyond this level, all iron assumed to be
% retained and removed from the system.

grafe = graze([A.idx.ps A.idx.pl],:,:) + pred([A.idx.ps A.idx.pl],:,:);

egefe = bsxfun(@times, grafe, [A.gs' 0 0]);
egefe = bsxfun(@times, egefe, permute(fe2n(:,[A.idx.ps A.idx.pl]), [2 3 1]));
egefe = sum(egefe, 2);

excfe = bsxfun(@times, grafe, [(1 - A.gs' - A.ge') 0 0]);
excfe = bsxfun(@times, excfe, permute(fe2n(:,[A.idx.ps A.idx.pl]), [2 3 1]));
excfe = sum(excfe,2);

lostfe = bsxfun(@times, grafe, [A.ge' 0 0]);
lostfe = bsxfun(@times, lostfe, permute(fe2n(:,[A.idx.ps A.idx.pl]), [2 3 1]));
lostfe = sum(lostfe,2);

egest(  [A.idx.psfe A.idx.plfe], A.idx.pofe, :) = egefe;
excrete([A.idx.psfe A.idx.plfe], A.idx.mys, :) = excfe;
graze(  [A.idx.psfe A.idx.plfe], A.idx.mys, :) = lostfe;