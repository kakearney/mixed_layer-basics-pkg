function [New, ped, old2new] = calcclusteredparams(Full, pedigree, idx, pedtype)
%CALCCLUSTEREDPARAMS Calculates new Ewe input params for clustered model
%
% [New, ped, old2new] = calcclusteredparams(Full, pedigree, idx, pedtype)
%
% Given an Ewe input structure and a list of clustering indices, this
% function calculates the new parameter values by summing biomass in each
% clustered group, and performing a biomass-weighted average on all other
% parameters.
%
% Input variables:
%
%   Full:       Ewe input structure for original food web
%
%   pedigree:   ng x 4 array pedigree values for B, P/B, Q/B, and DC
%
%   idx:        ng x 1 array of indices, where each index represents the
%               cluster number into which this group will be added.  Values
%               must be integers from 1 to nc, where nc is the number of
%               clusters,
%
%   pedtype:    method to use when calculating new pedigree
%               'propagate':    uses propagation of error through the
%                               summing/averaging calculations, assuming
%                               that pedigree represents 95% confidence
%                               interval (i.e. 2 std) in a normal
%                               distribution
%               'average':      performs a simple biomass-weighted average
%                               over clustered groups.  This will result in
%                               higher pedigree values than the propagation
%                               method, which may be preferable since the
%                               values are not based on measurements but
%                               rather rough guesses.
%
% Output variables:
%
%   New:        new Ewe input structure for clustered model, with nc groups
%
%   ped:        pedigree matrix for clustered model
%       
%   old2new:    ng x 1 array of indices showing which group each old one is
%               now part of (doesn't match idx because groups are resorted
%               by trophic level)
     


% Setup

[blah, names] = aggregate(idx, Full.name);
nfullcritter = cellfun('length', names);
hasmult = nfullcritter > 1;

EpFull = ecopathlite(Full);

%----------------------------
% Cluster pedigree-related 
% parameters (B, P/B, Q/B, 
% and DC) propagating error
%----------------------------

data = [EpFull.b EpFull.pb EpFull.qb];
stdev = (data .* pedigree(:,1:3))./2; % assume normal, 95% ci +/- percent = 2*std

% Biomass: sum biomass and propagate uncertainty:
% u{A+B} = sqrt(u{A}^2+u{B}^2)

[blah, bagg]    = aggregate(idx, data(:,1));
[blah, baggstd] = aggregate(idx, stdev(:,1));

bsum = cellfun(@sum, bagg);
bsumstd = cellfun(@(x) sqrt(sum(x.^2)), baggstd);

weights = cellfun(@(x) x./(sum(x)), bagg, 'uni', 0);

% P/B and Q/B weighted means and weighted standard deviation
% std(xavg) = sqrt(sum(std_i^2 .* w_i^2))

[blah, pbagg]    = aggregate(idx, data(:,2));
[blah, pbaggstd] = aggregate(idx, stdev(:,2));

pbavg = cellfun(@(x,y) sum(x.*y), pbagg, weights);
pbavgstd = cellfun(@(x,y) sqrt(sum(y.^2.*x.^2)), pbaggstd, weights);

[blah, qbagg]    = aggregate(idx, data(:,3));
[blah, qbaggstd] = aggregate(idx, stdev(:,3));

qbavg = cellfun(@(x,y) sum(x.*y), qbagg, weights);
qbavgstd = cellfun(@(x,y) sqrt(sum(y.^2.*x.^2)), qbaggstd, weights);

% Diet composition: Only given one value for each predator in the pedigree,
% so assume this fraction applies equally to all diet fractions

[blah, dcagg] = aggregate(idx, Full.dc');
[blah, dcaggstd] = aggregate(idx, pedigree(:,4));

nsimp = length(bsum);
dcavg = zeros(nsimp);
dcavgstd = zeros(nsimp);

for ig = 1:length(dcagg)
    dc1 = dcagg{ig}';
    ci1 = repmat(dcaggstd{ig}', size(dc1,1), 1);
    ci1(dc1 == 0) = 0;
    std1 = (dc1 .* ci1)./2;
    
    [blah, dc2] = consolidator(idx, dc1, @sum);
    [blah, std2] = consolidator(idx, std1, @(x) sqrt(sum(x.^2)));
    
    dc3 = sum(bsxfun(@times, dc2, weights{ig}'), 2);
    std3 = sqrt(sum(bsxfun(@times, std2.^2, weights{ig}'.^2), 2));
    
    dcavg(:,ig) = dc3;
    dcavgstd(:,ig) = std3;
   
end

%----------------------------
% New pedigree: propagated
%----------------------------

% For biomass, P/B, and Q/B, use the std values calculated via sum or
% weighted average

simpdata = [bsum pbavg qbavg];
simpstd = [bsumstd pbavgstd qbavgstd];

simppedi = (simpstd.*2)./simpdata;
is0 = simpstd == 0;
simppedi(is0) = 0;  % To avoid NaN

% For diet, can only have one value per predator.  

dcpedi = (dcavgstd.*2)./dcavg;
dcpedi = nanmean(dcpedi,1)';
dcpedi(isnan(dcpedi)) = 0.3; % For things that don't eat anything

simppedi = [simppedi dcpedi];

%----------------------------
% New pedigree: simple 
% average
%----------------------------

[blah, pedagg] = aggregate(idx, pedigree);
pedavg = cellfun(@(x,y) sum(bsxfun(@times,x,y),1), pedagg, weights, 'uni', 0);
pedavg = cell2mat(pedavg);

switch pedtype
    case 'propagated'
        ped = simppedi;
    case 'average'
        ped = pedavg;
end

%----------------------------
% Cluster all other 
% parameters using weighted 
% average method
%----------------------------

[blah, eeagg] = aggregate(idx, EpFull.ee);
eeavg = cellfun(@(x,y) sum(x.*y), eeagg, weights);

[blah, geagg] = aggregate(idx, EpFull.ge);
geavg = cellfun(@(x,y) sum(x.*y), geagg, weights);

[blah, baagg] = aggregate(idx, EpFull.ba);
baavg = cellfun(@(x,y) sum(x.*y), baagg, weights);

[blah, gsagg] = aggregate(idx, Full.gs);
gsavg = cellfun(@(x,y) sum(x.*y), gsagg, weights);

[blah, ppagg] = aggregate(idx, Full.pp);
ppavg = cellfun(@(x,y) sum(x.*y), ppagg, weights);

if ~all(floor(ppavg) == ppavg)
    warning('Not all pp values are integers');
end

ndet = sum(ppavg == 2);
if ndet > 1
    warning('More than one detritus group; you will have to manually adjust the detritus fate and discard fate parameters');
    multdet = true;
else
    multdet = false;
end

% Throw away "missing" values.  If all of the groups included in a cluster
% were missing a parameter, then I consider it missing. (Note that in the
% error propagation calcs, I used post-ecopath values to avoid any missing
% terms).

fun = @nanmean;
[blah, temp] = consolidator(idx, Full.b, fun);
nob  = isnan(temp);
[blah, temp] = consolidator(idx, Full.qb, fun);
noqb = isnan(temp);
[blah, temp] = consolidator(idx, Full.ee, fun);
noee = isnan(temp);
[blah, temp] = consolidator(idx, Full.ge, fun);
noge = isnan(temp);
[blah, temp] = consolidator(idx, Full.pb, fun);
nopb = isnan(temp);

bsum(nob)   = NaN;
qbavg(noqb) = NaN;
eeavg(noee) = NaN;
geavg(noge) = NaN;
pbavg(nopb) = NaN;

% Names

[blah, names] = aggregate(idx, Full.name);
names = cellfun(@(x) sprintf('%s\n', x{:}), names, 'uni', 0);
names = cellfun(@(x) x(1:end-1), names, 'uni', 0);

%----------------------------
% Build final food web
%----------------------------

ncluster = max(idx);

New.ngroup      = ncluster;
New.nlive       = ncluster - ndet;
New.ngear       = 1;
New.areafrac    = ones(ncluster,1);
New.b           = bsum;
New.pb          = pbavg;
New.qb          = qbavg;
New.ee          = eeavg;
New.ge          = geavg;
New.ba          = baavg;
New.gs          = gsavg;
New.dtImp       = zeros(ncluster,1);
New.bh          = New.b .* New.areafrac;
New.pp          = ppavg;
New.dc          = dcavg;
New.df          = ones(ncluster, ndet)./ndet;
New.immig       = zeros(ncluster,1);
New.emig        = zeros(ncluster,1);
New.emigRate    = zeros(ncluster,1);
New.ba          = zeros(ncluster,1);
New.baRate      = zeros(ncluster,1);
New.landing     = zeros(ncluster,1);
New.discard     = zeros(ncluster,1);
New.discardFate = ones(1, ndet)./ndet;            
New.maxrelpb    = ones(ncluster,1) * 2;
New.qbmaxqb0    = ones(ncluster,1) * 1000;
New.kv          = zeros(ncluster);
New.kv(New.dc>0) = 2;
New.name        = names;

% Reorder clustered groups by trophic level and primary production level

[New, srtidx] = sortewein(New);
ped = ped(srtidx,:);

[tf, old2new] = ismember(idx, srtidx);

