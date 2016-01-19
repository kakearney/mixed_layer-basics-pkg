function [S, A, Ep] = ecopathstats(A, varargin)
%ECOPATHSTATS Ecological network indices for Ecopath models
%
% S = ecopathstats(In)
% S = ecopathstats(In, p1, v1, ...)
%
% This function calculates a variety of network indices for an Ecopath
% model.  The indices calculated are based on the NetIndices R package;
% these are similar to (but in some cases not identical to) the formulas
% used by Statistics and Network Analysis tools in EwE6.
%
% Note: a few of these are still a work in progress, as noted in
% descriptions of output.  Don't use those values yet.  
%
% Input variables:
%
%   In:     ecopath input structure
%
% Optional input variables:
%
%   fleet:  How are fishing fleets treated?
%           'in':   Considered nodes of the system, catches are
%                   within-system fluxes, landings are considered exports
%                   from the fleet node and discards to detritus are
%                   within-system fluxes
%           'out':  Fleets are external to the system, the landing portion
%                   of catches is an export and the discard portion is
%                   routed directly from fish group to detritus 
%
%   x:      ensemble-related, see createensemble.m.  If included,
%           statistics will be calculated for each ensemble member where
%           approprate (indices based only on food web structure and not
%           parameter values, such as number of links, will be the same for
%           all and therefore will only be calculated once).
%
%   idx:    ensemble-related, see createensemble.m.  Indices corresonding
%           to x input.
%
% Output variables:
%
%   S:      n x 1 structure with the following fields.
%
%           T:          ng x ng array, Flow from compartment j to i, where
%                       j represents the columns of the flow matrix and i
%                       the rows (note that this row/column convention is
%                       the opposite of most of my Ecopath-related
%                       functions, but I've preserved it for consistency
%                       with the R package).
%
%           TSTp:       Total system throughput.  
%
%           TSTf:       Total system throughflow.
%
%           Ltot:       Number of links
%
%           Lint:       Number of internal links
%           
%           LD:         Link density
%
%           C:          Connectance
%
%           Tavg:       Average link weight
%
%           TSTfavg:    Average compartment throughflow
%
%           Cavg:       Compartmentalization
%
%           A:          Ascendency
%
%           DC:         Development capacity
%
%           O:          System overhead
%
%           AC:         Extent of development
%
%           TSTc:       Total system cycled throughflow (work in progress)
%
%           TSTs:       Total system non-cycled throughflow (work in
%                       progress)
%
%           FCI:        Finn's cycling index (work in progress)
%
%           FCIb:       Revised Finn's cycling index (work in progress)
%
%           Qsum:       Sum of consumption
%
%           Psum:       Sum of production
%
%           catchTL:    mean trophic level of catch
%
%           GE:         gross efficiency (catches/net primary production)

% Copyright 2016 Kelly Kearney

Opt.fleet = 'in';
Opt.x = [];
Opt.idx = [];

Opt = parsepv(Opt, varargin);

if isempty(Opt.x)
    Ep = ecopathlite(A, 'silent', true);
else
    [~, Ep] = ecopathlite(A, 'silent', true, 'x', Opt.x, 'idx', Opt.idx);
    A = subecopathens(A, Opt.x, Opt.idx);
end

% 
% if nargin < 2 || isempty(Ep)
%     Ep = ecopathlite(A, 'silent', true);
% end
% 

for ii = 1:length(Ep)

    %------------------
    % Transport/flow 
    % matrix
    %------------------

    switch Opt.fleet
        case 'in'
            if ii == 1
                exidx = [Ep(ii).Idx.out Ep(ii).Idx.lan];
                inidx = [Ep(ii).Idx.gpp];
                uuidx = Ep(ii).Idx.res;

                idx = setdiff(1:size(Ep(ii).flow,1), [exidx inidx uuidx]);
            end

            S(ii).n = length(idx);

            S(ii).T = zeros(S(ii).n + 3);
            S(ii).T(1:S(ii).n,1:S(ii).n) = Ep(ii).flow(idx,idx);

            S(ii).T(1:S(ii).n,S(ii).n+1) = sum(Ep(ii).flow(idx,exidx), 2);
            S(ii).T(1:S(ii).n,S(ii).n+2) = Ep(ii).flow(idx,uuidx);
            S(ii).T(S(ii).n+3,1:S(ii).n) = Ep(ii).flow(inidx,idx);

        case 'out'

            if ii == 1
                exidx = [Ep.Idx.out Ep.Idx.lan];
                inidx = [Ep.Idx.gpp];
                uuidx = Ep.Idx.res;

                idx = setdiff(1:size(Ep.flow,1), [exidx inidx uuidx Ep.Idx.gear]);
            end
            
            % First, reroute flows so bypass gears

            catches = [sum(A(ii).landing,2) A(ii).discard]; 
            fate1 = [zeros(1,A(ii).ngroup-A(ii).nlive); A(ii).discardFate]; 
            fate1 = [fate1 1-sum(fate1,2)];

            ttmp = Ep(ii).flow;
            ttmp(Ep(ii).Idx.gear,:) = 0;
            ttmp(:,Ep(ii).Idx.gear) = 0;
            ttmp([Ep(ii).Idx.liv Ep(ii).Idx.det], [Ep(ii).Idx.det, Ep(ii).Idx.lan]) = ...
                ttmp([Ep(ii).Idx.liv Ep(ii).Idx.det], [Ep(ii).Idx.det, Ep(ii).Idx.lan]) + catches*fate1;

            % Convert to T

            S(ii).n = length(idx);

            S(ii).T = zeros(S(ii).n + 3);
            S(ii).T(1:S(ii).n,1:S(ii).n) = ttmp(idx,idx);

            S(ii).T(1:S(ii).n,S(ii).n+1) = sum(ttmp(idx,exidx), 2);
            S(ii).T(1:S(ii).n,S(ii).n+2) = ttmp(idx,uuidx);
            S(ii).T(S(ii).n+3,1:S(ii).n) = ttmp(inidx,idx);

    end

    S(ii).T = S(ii).T'; % T_ij, i = sink, j = source

    %------------------------
    % General network indices
    %------------------------

    % Total system throughput

    S(ii).TSTp = sum(S(ii).T(:));

    % Total system throughflow

    S(ii).TSTf = sum(reshape(S(ii).T(1:S(ii).n,1:S(ii).n),[],1)) + ...
             sum(S(ii).T(:,S(ii).n+3)) - sum(S(ii).T(S(ii).n+1,:)) - sum(S(ii).T(S(ii).n+2,:));

    % Number of links

    if ii == 1
        S(ii).Ltot = nnz(S(ii).T);
    else
        S(ii).Ltot = S(1).Ltot;
    end

    % Number of internal links

    if ii == 1
        S(ii).Lint = nnz(S(ii).T(1:S(ii).n,1:S(ii).n));
    else
        S(ii).Lint = S(1).Lint;
    end

    % Link density

    if ii == 1
        S(ii).LD = S(ii).Ltot./S(ii).n;
    else
        S(ii).LD = S(1).LD;
    end

    % Connectance

    if ii == 1
        S(ii).C = S(ii).Lint./(S(ii).n*(S(ii).n-1));
    else
        S(ii).C = S(1).C;
    end

    % Average link weight

    S(ii).Tavg = S(ii).TSTp./S(ii).Ltot;

    % Average compartment throughflow

    S(ii).TSTfavg = S(ii).TSTf./S(ii).n; 

    % Compartmentalization

    if ii == 1
        adj = S(ii).T(1:S(ii).n,1:S(ii).n) > 0;
        c = zeros(size(adj));
        for ir = 1:S(ii).n
            for jj = 1:S(ii).n
                c(ir,jj) = (sum(adj(:,ir) & adj(:,jj)) + sum(adj(ir,:) & adj(jj,:))) ./ ...
                           (sum(adj(:,ir) | adj(:,jj)) + sum(adj(ir,:) | adj(jj,:)));
            end
        end

        S(ii).Cavg = sum(c(~eye(size(c))))./(S(ii).n*(S(ii).n-1));
    else
        S(ii).Cavg = S(1).Cavg;
    end

    %------------------------
    % Growth and development
    % Indices
    %------------------------

    % Ascendency

    TixTj = bsxfun(@times, sum(S(ii).T,1), sum(S(ii).T,2));

    tmp = S(ii).T .* log2(S(ii).T.*S(ii).TSTp./TixTj);
    S(ii).A = nansum(tmp(:));

    % Development capacity

    tmp = S(ii).T .* log2(S(ii).T./S(ii).TSTp);
    S(ii).DC = -nansum(tmp(:));

    % Overhead

    S(ii).O = S(ii).DC - S(ii).A;

    % Extent of development

    S(ii).AC = S(ii).A./S(ii).DC;

    %------------------------
    % Pathway analysis
    %------------------------

    % Total system cycled throughflow

    Tstar = S(ii).T(1:S(ii).n,1:S(ii).n);
    Tj = sum(Tstar,1); % sum outflows
    Ti = sum(Tstar,2); % sum inflows

    I = eye(S(ii).n);

    Gstar = Tstar./bsxfun(@max, Ti, Tj);
    Gstar(isnan(Gstar)) = 0;

    Q = inv(I - Gstar);

    S(ii).TSTc = sum((1 - diag(Q)).*Tj'); % TODO: not right, getting negatives
    
    % Total system non-cycled throughflow
    
    S(ii).TSTs = S(ii).TSTf - S(ii).TSTc;
    
    % Finn's cycling index
    
    S(ii).FCI = S(ii).TSTc./S(ii).TSTf;
    
    % Revised Finn's cycling index
    
    S(ii).FCIb = S(ii).TSTc./S(ii).TSTp;

    %------------------------
    % Environ analysis
    %------------------------

%     % Transitive closure matrix
%     
%     S(ii).G = bsxfun(@rdivide, Tstar, Tj);
%     S(ii).G(isnan(S(ii).G)) = 0;
%     
%     % Integral nondimensional matrix
%     
%     S(ii).N = inv(I - S(ii).G);

    %-------------------
    % Ecopath statistics
    %-------------------
    
    S(ii).Qsum = sum(Ep(ii).q0sum(1:A(ii).nlive)); % Sum of consumption
    S(ii).Psum = sum(Ep(ii).pb .* Ep(ii).b); % Sum of production
    
    fmort = Ep(ii).fishMortRate .* Ep(ii).b;
    S(ii).catchTL = sum(fmort./sum(fmort) .* Ep(ii).trophic);
    
    netpp = sum(Ep(ii).pb .* Ep(ii).b .* (A(ii).pp == 1));
    S(ii).GE = sum(fmort)./netpp;
    

end




