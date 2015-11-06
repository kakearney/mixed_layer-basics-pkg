function Ens = createensemble(Ewein, pedigree, nsample, varargin)
%CREATEENSEMBLE Generate ensemble of Ecopath models
%
% Set = createensemble(Ewein, pedigree, nset)
%
% This function creates an ensemble of Ecopath models.  It was designed as
% a method to include uncertainty in Ecopath-derived models (namely, my
% water column ecosystem model, which couples an Ecosim-like food web model
% to a biogeochemical model) through the generation of a large number of
% potential Ecopath models that all fall within the prescribed uncertainty
% ranges, rather than relying on a single Ecopath parameter set.
%
% It may also be a useful tool for building a balanced model.  This
% function does not require that the original Ecopath model included
% (Ewein) itself be balanced.  Given an unbalanced model along with a
% pedigree array, this function will be able to find any parameter sets
% within those confidence bounds that do balance.
%
% A brief discussion of the use of this method to parameterize a model can
% be found here:
%
% Kearney, K. A., Stock, C., Aydin, K. & Sarmiento, J. L. Coupling
% planktonic ecosystem and fisheries food web models for a pelagic
% ecosystem: Description and validation for the subarctic Pacific.
% Ecological Modelling 237-238, 43?62 (2012).  
%
% For a more in-depth discussion of the process, see Chapter 2 in
%
% "An analysis of marine ecosystem dynamics through development of a
% coupled physical-biogeochemical-fisheries food web model", Ph.D.
% dissertation, Kelly Kearney, Princeton University
%
% Update: 04/02/2014...
%
% I've now added the ability to use latin hypercube sampling, rather than
% Monte Carlo sampling, to generate the point sets. I found that this
% sampling technique presents a way to maximize the amount of parameter
% space covered by a relatively small number of samples in a very
% high-dimensional parameter space.  I experimented with a few other
% techniques, such as centroidal Voronoi tessellation, which in theory
% should produce a more uniform distribution of points, but this resulted
% in poor filling of the univariate parameter space simply due to the
% typically small number of samples relative to the very high-dimensional
% parameter space... plus it gets numerically tricky to maintain the
% uniformity once I add in the mass-balance criterion. So for now, I've
% only included the two sampling techniques.
%
% Input variables:
%
%   Ewein:      Ewe input structure (see ecopathlite.m for details)
%
%   pedigree:   ngroup x 6 array of pedigree values (B, PB, QB, DC, EE,
%               GE).  All values should be between 0 and 1.  The last two
%               columns (EE, GE) are optional; if not included, they will
%               be assigned NaNs.  The pedigree values assign uncertainty
%               to the parameter values (see pdfname input for specifics).
%               A value of 0 or NaN indicates no variation, and
%               corresponding parameters will be left as point estimates.
%
%   nset:       number of ensemble members to collect
%
% Optional input variables (passed as parameter/value pairs):
%
%   collect:    'all' or 'balanced', specifies whether to return all
%               parameter sets, or only those that balance.  In the latter
%               case, parameter sets will be generated until nset balanced
%               sets are found; this may take a while, depending on the
%               particulars of the input ecosystem. ['balanced']
%
%   pdfname:    Specifies the probability distribution from which
%               parameters will be chosen. Note that final diet composition
%               will always be renormalized such that predator diet sums to
%               1, resulting in a normal distribution for each diet portion
%               across ensemble members, but the pre-normalized choices
%               will follow the specified distribution. (Note: I used to
%               offer a truncated Gaussian option, but that option isn't a
%               "real" distribution, making it difficult to keep that part
%               of the code in sync with the rest... that option was mostly
%               for some parameter tests anyway, so I've now eliminated it)
%               ['uniform']
%
%               'uniform':  Uniform distribution between x-ped*x and
%                           x+ped*x, where x is the point value from the
%                           original Ecopath model and ped is the pedigree
%                           value.
%               'lognormal':Lognormal distribution with mean x and variance
%                           (ped*x/2)^2.
%
%   sample:     Specifies the technique used to choose samples from the
%               parameter distribution functions. ['mcs']
%
%               'mcs':      Monte Carlo simulation
%
%               'lhs':      Latin hypercube simulation, with a maximin
%                           criterion to maximize point distance between
%                           ensemble members
%
%   lhsiter:    For latin hypercube sampling only, number of iterations to
%               perform per sample block in an attempt to minimize distance
%               between points (i.e. 'iterations' parameter to pass to
%               lhsdeign.m) [20]
%
% Output variables:
%
%   Set:        1 x 1 structure array with the following fields:
%
%               Ewein:  the Ewe input stucture used to generate this
%                       ensemble
%
%               idx:    1 x 6 cell array, corresponding to biomass (B),
%                       production/biomass (PB), consumption/biomass (QB),
%                       diet composition (DC), ecotrophic efficiency (EE),
%                       and growth efficiency (GE).  The cell arrays hold
%                       the indices of the groups (or in the case of DC,
%                       the linear array index of the prey/predator links)
%                       for which ensemble values were generated (i.e.
%                       non-NaN or non-zero values, as appropriate to the
%                       given parameter, in the original model).
%                      
%               x:      1 x 6 cell array of parameter values.  The arrays
%                       in each cell are nset x nidx arrays, where nidx
%                       corresponds to the length of the cells in idx.
%                       Note that these values include the
%                       randomly-distributed values *without* any
%                       inter-parameter adjustments (such as DC
%                       normalization or multi-stanza growth curve fits);
%                       these adjustments must be applied before being used
%                       in a new Ecopath model (both subecopathens.m and
%                       ecopathlite.m with ensemble input option perform
%                       these adjustments).
%
%               mu:     1 x 6 cell array of mu values associated with the
%                       lognormal distribution for each ensemble set
%
%               sig:    1 x 6 cell array of sigma values associated with
%                       the lognormal distribution of each ensemble set
%
%               lo:     1 x 6 cell array of the lower-end cutoff for the
%                       uniform distribution of each ensemble set
%
%               hi:     1 x 6 cell array of the upper-end cutoff for the
%                       uniform distribution of each ensemble set
%
%               nall:   scalar, number of total parameter sets generated in
%                       order to create the required number of parameter
%                       sets.  For collect = 'all', this will be the same
%                       as nset; for collect = 'balanced', it may be much
%                       higher.

% Copyright 2012-2014 Kelly Kearney
% kakearney@gmail.com


% Parse optional inupts

Opt.pdfname = 'uniform';
Opt.collect = 'balanced';
Opt.maxiter = Inf;
Opt.sample  = 'mcs';
Opt.lhsiter = 20;

Opt = parsepv(Opt, varargin);

% Number of groups shortcut

ng = Ewein.ngroup;

%----------------
% Gather vars
%----------------

pedigree(pedigree == 0) = NaN;

bbped = pedigree(:,1);
pbped = pedigree(:,2);
qbped = pedigree(:,3);
dcped = repmat(pedigree(:,4), 1, ng)';
if size(pedigree,2) > 4
    eeped = pedigree(:,5);
else
    eeped = nan(size(bbped)); % For back-compatability, more or less
end
if size(pedigree,2) > 5
    geped = pedigree(:,6);
else
    geped = nan(size(bbped));
end

bbmid = Ewein.b;
pbmid = Ewein.pb;
qbmid = Ewein.qb;
dcmid = Ewein.dc;
eemid = Ewein.ee;
gemid = Ewein.ge;

bbidx = find(~isnan(bbmid) & ~isnan(bbped));
pbidx = find(~isnan(pbmid) & pbmid ~= 0 & ~isnan(pbped));
qbidx = find(~isnan(qbmid) & qbmid ~= 0 & ~isnan(qbped));
dcidx = find(dcmid & ~isnan(dcped));
eeidx = find(~isnan(eemid) & Ewein.pp ~= 2 & ~isnan(eeped));
geidx = find(~isnan(gemid) & any(isnan([pbmid qbmid]),2) & ~isnan(geped)); % GE overwritten if PB and QB also provided

idx = {bbidx, pbidx, qbidx, dcidx, eeidx, geidx};

% All uncertain variables, together

vmid = [bbmid(bbidx); pbmid(pbidx); qbmid(qbidx); dcmid(dcidx); eemid(eeidx); gemid(geidx)];
vped = [bbped(bbidx); pbped(pbidx); qbped(qbidx); dcped(dcidx); eeped(eeidx); geped(geidx)];
vciv = vmid .* vped;

nvar = length(vmid);

% Convert to mu/sigma values (for lognormal) ...

vvar = (vciv./2).^2; % variance

mu = log((vmid.^2)./sqrt(vvar+vmid.^2));
sigma = sqrt(log(vvar./(vmid.^2)+1));

% ... or upper/lower (for uniform)

lo = vmid - vciv;
hi = vmid + vciv;

%----------------
% Generate sets
%----------------

S = warning('off', 'EWE:cannibalTooHigh');

switch Opt.collect
    
    case 'all'
        
        % Sample points uniformly in CDF
        
        switch Opt.sample
            case 'mcs'
               p = rand(nsample, nvar);
            case 'lhs'
               p = lhsdesign(nsample, nvar, ...
                   'criterion', 'maximin', ...
                   'iterations', Opt.lhsiter);
            otherwise
                error('Unrecognized sample option');
        end
       
        % Translate back to values based on each CDF
    
        switch Opt.pdfname
            
            case 'uniform'
                x = zeros(size(p));
                for iv = 1:size(p,2)
                    x(:,iv) = unifinv(p(:,iv), lo(iv), hi(iv));
                end
                
            case 'lognormal'
                x = zeros(size(p));
                for iv = 1:size(p,2)
                    x(:,iv) = logninv(p(:,iv), mu(iv), sigma(iv));
                end
                
            otherwise
                error('Unrecognized distribution');
        end
        nall = nsample;
        
    case 'balanced'
        
        nbal = 0;
        nall = 0;
        x = zeros(nsample, nvar);
        
        
        cpb = ConsoleProgressBar();
        cpb.setMaximum(nsample);
        fprintf('Generating ensembles...\n');
        cpb.start();

        while nbal < nsample
            
            % Sample points uniformly in CDF

            switch Opt.sample
                case 'mcs'
                   ptmp = rand(nsample, nvar);
                case 'lhs'
                   ptmp = lhsdesign(nsample, nvar, ...
                       'criterion', 'maximin', ...
                       'iterations', 20);
                otherwise
                    error('Unrecognized sample option');
            end

            % Translate back to values based on each CDF

            switch Opt.pdfname

                case 'uniform'
                    xtmp = zeros(size(ptmp));
                    for iv = 1:size(ptmp,2)
                        xtmp(:,iv) = unifinv(ptmp(:,iv), lo(iv), hi(iv));
                    end

                case 'lognormal'
                    xtmp = zeros(size(ptmp));
                    for iv = 1:size(ptmp,2)
                        xtmp(:,iv) = logninv(ptmp(:,iv), mu(iv), sigma(iv));
                    end

                otherwise
                    error('Unrecognized distribution');
            end

            % Keep only balanced sets
            
            [isbal, xtmp] = checkbalance(Ewein, idx, xtmp);

            xbal = xtmp(isbal,:);
            nnew = size(xbal,1);

            if nnew+nbal <= nsample
                x((1:nnew)+nbal,:) = xbal;
                nall = nall + nsample;
            else
                x((nbal+1):end,:) = xbal(1:(nsample-nbal),:);
                nall = nall + (nsample-nbal);
            end

            nbal = nbal + nnew;

            cpb.setValue(min(nbal,nsample));
            cpb.setText(sprintf('Attempted: %d', nall));
        end
        cpb.stop();
        fprintf('\n');        
end

%----------------
% Reformat for 
% output
%----------------

warning(S);

% Additional parameters

Ens.Ewein = Ewein;
Ens.idx   = {bbidx, pbidx, qbidx, dcidx, eeidx, geidx};
ncol = cellfun(@length, Ens.idx);
Ens.x     = mat2cell(x, nsample, ncol);
Ens.mu    = mat2cell(mu,    ncol);
Ens.sig   = mat2cell(sigma, ncol);
Ens.lo    = mat2cell(lo,    ncol);
Ens.hi    = mat2cell(hi,    ncol);
Ens.nall  = nall;

%*********************

%----------------
% Check for mass
% balance
%----------------

function [isbal,x] = checkbalance(Ewein, idx, x) %bbidx, pbidx, qbidx, dcidx, x)

nx = size(x,1);
ncol = cellfun(@length, idx);
x = mat2cell(x, nx, ncol);

[blah, Ep] = ecopathlite(Ewein, 'x', x, 'idx', idx, 'skipextra', true, 'silent', true);

isbal = all([Ep.ee] <= 1 & [Ep.ee] >= 0 & ~isnan([Ep.ee]) & ...
            [Ep.ge] >= 0 & [Ep.ge] <= 1, 1);

x = cat(2, x{:});

    
    

