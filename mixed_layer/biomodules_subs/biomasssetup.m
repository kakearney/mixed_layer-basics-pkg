function [bv, ba, basum, bfrac, zlfrac, nb, nz] =  biomasssetup(bio, A)
%BIOMASSSETUP Biomass calculations for nemurokak/wce modules
%
% [bv, ba, basum, bfrac, zlfrac, nb, nz] =  biomasssetup(bio, A)
%
% The diapause and preyvis options require different processes to "see"
% different fractions of the various functional group biomasses.
%
% Input variables:
%
%   bio:    nz x nbsv array of biomass, values received from ode solver
%
%   A:      structure various wce/nemurokak parameters
%
% Output variables:
%
%   bv:     volumetric, per-layer biomass concentration
%           orig:       nz x nb, same as input
%           zlcombo:    nz x nb, ZL1/ZL2 moved to ZL (diapause only)
%           prey:       nb x nb x nz. prey by pred by depth, prey profile
%                       for each predator/prey link
%           pred:       nb x nb x nz. prey by pred by depth, predator
%                       profile for each predator/prey link 
%
%   ba:     depth-integrated, per-layer biomass, same fields as bv
%
%   basum:  water column-integrated biomass, same fields as bv
%
%   zlfrac: fraction of ZL1 in the ZL total
%
%   nb:     number of biological state variables
%
%   nz:     number of depth layers

% Copyright 2014 Kelly Kearney

%---------------------
% Volumetric biomass
% (per layer, mol/m^3)
%---------------------

% Biomass as passed by mixed layer.  In the diapause case, ZL is 0 and ZL1
% and ZL2 hold the non-diapausing and diapausing copepod populations,
% respectively.

bv.orig = bio; 

% The only alteration here is to move the split ZL1/ZL2 populations into ZL
% if diapause is on, and zero out those subgroups.

bv.zlcombo = bv.orig;
if A.diapause
    bzl = bv.orig(:,[A.idx.zl1 A.idx.zl2]);
    
    bv.zlcombo(:,A.idx.zl) = sum(bzl,2);
    bv.zlcombo(:,[A.idx.zl1 A.idx.zl2]) = 0;
    
    zlfrac = bsxfun(@rdivide, bzl, sum(bzl,2));
    
end

% The predator/prey population varies by diet link, since some predators
% see prey differently.  ZL1/ZL2 is seen as a combined group, and predator
% ability to access prey is set per depth layer.

[nz, nb] = size(bv.orig);

[bv.pred, bv.prey] = deal(zeros(nb,nb,nz));
for iz = 1:nz
    btmp = bv.zlcombo(iz,:)';
    
    bv.prey(:,:,iz) = btmp * A.preyvis(iz,:);
    bv.pred(:,:,iz) = ones(nb,1) * btmp';
    
end

% TODO: make sure nekton see other nekton?
    

%---------------------
% Integrated biomass
% (per layer, mol/m^3)
%---------------------

[ba.orig,    basum.orig,    bfrac.orig]    = intoverdepth(bv.orig, A.dz, 1);
[ba.zlcombo, basum.zlcombo, bfrac.zlcombo] = intoverdepth(bv.orig, A.dz, 1);
[ba.pred,    basum.pred,    bfrac.pred]    = intoverdepth(bv.pred, A.dz, 3);
[ba.prey,    basum.prey,    bfrac.prey]    = intoverdepth(bv.prey, A.dz, 3);

%---------------------
% Subfunction: 
% integrate over depth
%---------------------   

function [ba, basum, bfrac] = intoverdepth(bv, dz, dim)
                     
if dim == 1
    dz = dz;
elseif dim == 3
    dz = permute(dz, [2 3 1]);
end
    
ba    = bsxfun(@times, bv, dz);       % nz x nbsv, mol N/m^2 per layer
basum = sum(ba, dim);                 % 1 x nbsv, mol N/m^2 total water column
bfrac = bsxfun(@rdivide, ba, basum);  % fraction of biomass in each layer

