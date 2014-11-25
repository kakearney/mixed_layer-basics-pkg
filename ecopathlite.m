function varargout = ecopathlite(S)
%ECOPATHLITE Rewrite of Ecopath algorithms
%
% ecopathlite(S)
% C = ecopathlite(S)
% [C, flag] = ecopathlite(S)
%
% This function reproduces the main calculations performed by the Ecopath
% portion of the EwE model (www.ecopath.org).  If no output variable is
% provided, the results are printed to tables in the command window. 
%
% Ecopath is used to calculate a snapshot of an ecosystem, including the
% biomass of all functional groups (living groups, fishing gears, and
% detrital pools) and the fluxes between these groups.  This function is a
% bare-bones version of the algorithm; it is not really meant to be used to
% set up and balance a model for the first time, since it does not provide
% any of the visual checks on the results (e.g. whether EE values are > 1,
% etc); use the original Ecopath software if you're looking for this type
% of behavior.  Nor does it include many of the more complicated setup
% options: multi-stanza groups, economic variables (e.g. market prices,
% fleet dynamics), etc.  It simply provides the initial mass-balance
% calculations, in a stripped-down, easy-to-automate format.
%
% The units below are defined in terms of three variables: M = mass, A =
% area (or volume), and T = time. In the original software, the default is
% M = metric tons wet weight, A = km^2, and T = year.  You can use whatever
% units work best for your own purposes, as long as they remain consistent
% across all variables.  
%
% Input and output variables were chosen to be somewhat consistent with the
% tables seen in Ecopath with Ecosim version 5.
%
% For more information on the Ecopath concept, see:
%
% Christensen, V. & Pauly, D. ECOPATH II--a software for balancing
% steady-state ecosystem models and calculating network characteristics.
% Ecological Modelling 61, 169?185 (1992).  
%
% Christensen, V. & Walters, C. J. Ecopath with Ecosim: methods,
% capabilities and limitations. Ecological Modelling 172, 109?139 (2004). 
%
% Note: I developed this code based on the equations documented in
% Appendix 4 of the EwE5 help files (this appendix is referenced in the EwE
% User's Guide, but only seems to be available through the help menu of
% EwE5).
%
% Input variables:
%
%   S:  structure with the following fields.  Values of the fields b, pb,
%       qb, ee, ge, gs, and/or dtImp that are defined as NaN indicate
%       unknown values, which will be filled in by the ecopath algorithm.
% 
%       ngroup:         1 x 1 array, number of functional groups in the
%                       model
%
%       nlive:          1 x 1 array, number of live (non-detrital) groups
%                       in the model  
%
%       ngear:          1 x 1, number of fishing gear types in the model
%
%       areafrac:       ngroup x 1 array, fraction of habitat area occupied
%                       by each group (no units, 0-1)
%
%       b:              ngroup x 1 array, biomass (M A^-1)
%
%       pb:             ngroup x 1 array, production over biomass ratios
%                       (T^-1) 
%   
%       qb:             ngroup x 1 array, consumption over biomass ratios
%                       (T^-1) 
%
%       ee:             ngroup x 1 array, ecotrophic efficiencies (no
%                       units, 0-1) 
% 
%       ge:             ngroup x 1 array, gross efficiency, i.e. production
%                       over consumption ratio (no units)
%
%       gs:             ngroup x 1 array, fraction of consumed food that is
%                       not assimilated (no units)
%
%       dtImp:          ngroup x 1 array, detritus import (should be zero
%                       for all non-detrital groups) (M A^-1 T^-1)
%
%       bh:             ngroup x 1 array,  habitat biomass, i.e. biomass
%                       per unit habitable area (M A^-1). This variable
%                       is designed as a shortcut if all your critters are
%                       clustered in only a small portion of the habitat
%                       area, such that bh = b/areafrac.
%
%       pp:             ngroup x 1 array, fraction of diet consisting of
%                       primary production, pp = 2 indicates detritus
% 
%       dc:             ngroup x ngroup array, diet composition, dc(i,j)
%                       tells fraction predator j's diet consisting of prey
%                       i (no units, 0-1)
% 
%       df:             ngroup x (ngroup - nlive) array, fraction of each
%                       group that goes to each detrital group due to other
%                       mortality and egestion (no units, 0-1)
%
%       immig:          ngroup x 1 array, immigration into area 
%                       (M A^-1 T^-1) 
%
%       emig:           ngroup x 1 array, emigration out of area 
%                       (M A^-1 T^-1)
%
%       emigRate:       ngroup x 1 array, emigration per unit biomass
%                       (T^-1) 
% 
%       ba:             ngroup x 1 array, biomass accumulation 
%                       (M A^-1 T^-1)
%
%       baRate          ngroup x 1 array, biomass accumulation per unit
%                       biomass (T^-1)
%
%       landing         ngroup x ngear array, landings of each group by
%                       each gear type (M A^-1 T^-1 ?)
%
%       discard         ngroup x ngear array, discards of each group by
%                       each gear type (M A^-1 T^-1 ?)
%
%       discardFate:    ngear x (ngroup - nlive) array, fraction of
%                       discards from each gear type that go to each
%                       detritus group (no units, 0-1)
%
% Output variables:
%
%   C:  structure with the following fields, all ngroup x 1 arrays unless
%       otherwise specified
%
%       trophic:        ngroup x 1 array, trophic level of each group (no
%                       unit) 
%
%       areafrac:       ngroup x 1 array, fraction of total area occupied
%                       by group (no units, 0-1)
%
%       bh:             ngroup x 1 array, biomass in habitable area 
%                       (M A^-1)
%
%       b:              ngroup x 1 array, total biomass (M A^-1)
%
%       pb:             ngroup x 1 array, production/biomass ratio (T^-1)
%
%       qb:             ngroup x 1 array, consumption/biomass ratio (T^-1)
%
%       ee:             ngroup x 1 array, ecotrophic efficiency (no unit)
%
%       ge:             ngroup x 1 array, growth efficiency, i.e.
%                       production/consumption (no unit)
%
%       ba:             ngroup x 1 array, biomass accumulation 
%                       (M A^-1 T^-1)
%
%       baRate:         ngroup x 1 array, biomass accumulation per unit
%                       biomass (T^-1) 
%
%       migration:      ngroup x 1 array, net migration (M A^-1 T^-1)
%
%       flowtodet:      (ngroup + ngear) x 1 array, flow to detritus from
%                       each group and each gear type (M A^-1 T^-1)
%
%       fishMortRate:   ngroup x 1 array, mortality per unit biomass due to
%                       fishing (T^-1)
%
%       predMortRate:   ngroup x 1 array, mortality per unit biomass due to
%                       predation, M2 in some documentation (T^-1) 
%
%       migrationRate:  ngroup x 1 array, net migration per unit biomass
%                       (T^-1)
%
%       otherMortRate:  ngroup x 1 array, mortality per unit biomass due to
%                       anything else (T^-1)
%
%       predMort:       ngroup x ngroup array, predation mortality broken
%                       down by predator and prey, C.predMoreRate =
%                       sum(C.predMort,2). (T^-1)
%
%       q0:             ngroup x ngroup array, q0(i,j) is the flux of
%                       biomass from group i to group j.  For j = live,
%                       this is due to predation.  For j = detrital, this
%                       includes non-predatory mortality and egestion. 
%                       (M A^-1 T^-1)
%
%       q0sum:          ngroup x 1 array, total consumption by each
%                       predator (M A^-1 T^-1)
%
%       respiration:    ngroup x 1 array, respiration (M A^-1 T^-1)
%
%       searchRate:     ngroup x ngroup array, search rates of each
%                       predator for each prey, assuming simple linear
%                       dynamics, i.e. Qij = a*Bi*Bj (A M^-1 T^-1) 
%
%       detexport:      (ngroup-nlive) x 1 array, amount of detritus
%                       exported from the system (M A^-1 T^-1)
%
%   flag:   false if all unknowns are filled, true if not.  This output was
%           primarily added for testing purposes.  Unfilled unknowns are
%           usually due to incorrect input, but may point to a bug in this
%           implementation of the Ecopath algorithm.

% Copyright 2012 Kelly Kearney
% kakearney@gmail.com

%------------------------------
% Setup
%------------------------------

if ~isscalar(S)
    error('Input structure must be scalar');
end

S = ecopathinputcheck(S, true);

%------------------------------
% Setup calculations
%------------------------------

islive = (1:S.ngroup)' <= S.nlive; % Logical mask for live groups

% If emigration and biomass is given as an input rather than emigration per
% unit biomass (i.e. emigration rate), calculate the rate 

temp = (S.emig(islive) > 0) & ~isnan(S.b(islive)) & (S.emigRate(islive) == 0);
S.emigRate(temp) = S.emig(temp) .* S.b(temp);

% Calculate growth efficiency (i.e. production/consumption ratio),
% production/biomass ratio, and consuption/biomass ratio if needed and
% possible

[S.pb, S.qb, S.ge] = pbq(S.pb, S.qb, S.ge, islive);

% Calculate total catches for each group

catches = sum(S.landing + S.discard, 2);

%--------------------------
% Determine some predator/
% prey relationships
%--------------------------

% Prey and predator masks

isprey = S.dc > 0;
ispred = S.dc' > 0;

% Predators not including group itself

isPredNoCannib = ispred & ~eye(S.ngroup);

%--------------------------
% Algorithms to calculate
% missing variables
%--------------------------

flag = false;
while ~checkbasic(S.b, S.pb, S.qb, S.ee, S.ge, islive, S.pp)

    param = [S.b S.pb S.qb S.ee S.ge];
    
    % Run the p-b-q algebra again (in case new values have been filled in)
    
    [S.pb, S.qb, S.ge] = pbq(S.pb, S.qb, S.ge, islive);
    
    % Algorithm 1: Estimation of P/B

    knowPredInfo = ~any(bsxfun(@and, ispred, isnan(S.b))) & ...
                   ~any(bsxfun(@and, ispred, isnan(S.qb)));
    
    m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2); 

    nm = calcnetmigration(S.emigRate, S.emig, S.immig);
    
    canRun = (islive & ...          % only live groups           
              isnan(S.pb) & ...     % P/B unknown
              ~isnan(S.b) & ...     % B known
              ~isnan(S.ee) & ...    % EE known
              knowPredInfo');       % know B and Q/B for all group's predators
      
    S.pb(canRun) = (catches(canRun) + nm(canRun) + S.ba(canRun) + ...
                    m2(canRun))./(S.b(canRun) .* S.ee(canRun)); 
    
    % Algorithm 2: Estimation of EE
    
    knowPredInfo = ~any(bsxfun(@and, ispred, isnan(S.b))) & ...
                   ~any(bsxfun(@and, ispred, isnan(S.qb)));
    
    m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2); 
    nm = calcnetmigration(S.emigRate, S.emig, S.immig);
    
    canRun = (islive & ...          % only live groups           
              isnan(S.ee) & ...     % EE unknown
              ~isnan(S.b) & ...     % B known
              ~isnan(S.pb) & ...    % P/B known
              knowPredInfo');       % know B and Q/B for all group's predators

    S.ee(canRun) = (catches(canRun) + nm(canRun) + S.ba(canRun) + ...
                    m2(canRun))./(S.b(canRun) .* S.pb(canRun)); 
    
    % Algorithm 3: Dealing with B and Q/B as unknowns (Note: never found
    % model that needs this, so never tested)
    
    knowPredInfo = ~any(bsxfun(@and, isPredNoCannib, isnan(S.b))) & ...
                   ~any(bsxfun(@and, isPredNoCannib, isnan(S.qb)));
    
               
    knowAllPreyInfo = bsxfun(@and, isprey, ~isnan(S.b)) & ...
                      bsxfun(@and, isprey, ~isnan(S.pb)) & ...
                      bsxfun(@and, isprey, ~isnan(S.ee)) & ...
                      bsxfun(@and, isprey, knowPredInfo);
    
    
    canRun = (islive & ...                          % only live groups
              (( isnan(S.b) & ~isnan(S.qb)) | ...   % B unknown & Q/B known
               (~isnan(S.b) &  isnan(S.qb))) & ...  % or B known & Q/B unknown
              knowPredInfo' & ...                   % know B and Q/B of predators except itself
              any(knowAllPreyInfo, 1)');     % know B, P/B, EE, and pred info (except group) for at least one prey
          
    group1 = canRun & isnan(S.b); % B unknown groups
    group2 = canRun & ~group1;    % Q/B unknown groups
    
    ng = length(S.b);
    k = zeros(ng,1);
    idxknowprey = find(any(knowAllPreyInfo, 1));
    for ii = idxknowprey
        k(ii) = find(knowAllPreyInfo(:,ii), 1, 'first');
    end
    iitmp = (1:ng)';
    idx = ones(ng,1); % This is just to avoid 0 subscripts in testing; we won't actually use the placeholders
    idx(idxknowprey) = sub2ind([ng ng], k(idxknowprey), iitmp(idxknowprey));
    
    partm2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc .* ~eye(S.ngroup)), 2); % predation, not including cannibalism      
    m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2); % all predation
    nm = calcnetmigration(S.emigRate, S.emig, S.immig);
    dcii = diag(S.dc);
    dcik = S.dc(idx);
    
    S.b(group1) = partm2(group1) + catches(group1) + nm(group1) + ...
                  dcii(group1) .*(S.b(k(group1)) .* S.pb(k(group1)) .* ...
                  S.b(k(group1)) .* S.ee(k(group1)) - catches(k(group1)) - ...
                  nm(k(group1)) - m2(k(group1)))./dcik(group1);
    
    S.qb(group2) = (S.b(k(group2)) .* S.pb(k(group2)) .* S.b(k(group2)) .* ...
                   S.ee(k(group2)) - catches(k(group2))) ./ ...
                   (dcik(group2) ./ S.b(group2));
    
    % Algorithm 4: Estimating biomasses only
    
    knowPredInfo = ~any(bsxfun(@and, isPredNoCannib, isnan(S.b))) & ...
                   ~any(bsxfun(@and, isPredNoCannib, isnan(S.qb)));
    
    
    canRun = (islive & ...          % only live groups
              isnan(S.b) & ...      % B unknown
              ~isnan(S.pb) & ...    % P/B known
              ~isnan(S.ee) & ...    % EE known
              ~isnan(S.qb) & ...    % Q/B known
              knowPredInfo');       % B and Q/B of predators except itself known
    
    partm2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc .* ~eye(S.ngroup)), 2); % predation, not including cannibalism      
    nm = calcnetmigration(S.emigRate, S.emig, S.immig);
    dcCannib = diag(S.dc);
    
    cannibCheck1 = canRun & (S.pb .* S.ee) == (S.qb .* dcCannib);
    cannibCheck2 = canRun & (S.pb .* S.ee) < (S.qb .* dcCannib);
    if any(cannibCheck1)
        idx = find(cannibCheck1);
        idxStr = sprintf('%d,', idx);
        idxStr = idxStr(1:end-1);
        warning('EWE:cannibalWithoutB', 'Group(s) (%s) are missing biomass but are only preyed on by themselves; group(s) must be split in two to solve\nExiting without solution', idxStr);
        break
    end
    if any(cannibCheck2)
        idx = find(cannibCheck2);
        idxStr = sprintf('%d,', idx);
        idxStr = idxStr(1:end-1);
        warning('EWE:cannibalTooHigh', 'Group(s) (%s) have cannibalism losses that exceed predation mortality\nExiting without solution', idxStr); 
        break
    end
    
    S.b(canRun) = (catches(canRun) + nm(canRun) + S.ba(canRun) + ...
                   partm2(canRun)) ./ (S.pb(canRun) .* S.ee(canRun) - ...
                   S.qb(canRun) .* dcCannib(canRun));
               
    % Algorithm 5: The generalized inverse
    
    % TODO: Not really clear on when this is used, or how.  Is X the same
    % as B in these equations?  And is this a necessary step, or just a
    % rephrasing of algorithms 3 & 4?
    
    % If we made it through all algorithms without filling in any new
    % numbers, break out of loop
    
    if isequalwithequalnans(param, [S.b S.pb S.qb S.ee S.ge])
        warning('EWE:unknownsremain', 'Unable to fill in all unknowns');
        flag = true;
        break
    end
    
end

needBh = isnan(S.bh) & ~isnan(S.b);
S.bh(needBh) = S.b(needBh) ./ S.areafrac(needBh);

%--------------------------
% Set some values to 0 that
% haven't been corrected
%--------------------------

S.pb(S.pp == 2) = 0;    % No production for detritus
S.qb(S.pp >= 1) = 0;    % No consumption for primary producers
S.ge(S.pp >= 1) = 0;    % Q = 0 for primary producers, so P/Q = Inf, set to 0 instead just as placeholder
S.gs(S.pp >= 1) = 0;    % Q = 0 so unassim irrelevant for detritus and producers

%--------------------------
% Detritus calculations
%--------------------------

% Detritus produced from mortality and egestion

mort = bsxfun(@times, S.b .* S.pb .* (1 - S.ee), S.df);   
egest = bsxfun(@times, S.b .* S.qb .* S.gs, S.df);   
detgroups = mort + egest;

% Detritus produced from fisheries discards

detfisheries = bsxfun(@times, sum(S.discard,1)', S.discardFate);

% Input to detritus groups from groups, fleets, and import

det = [detgroups; detfisheries];            % By source and destination
flowtodet = sum(det, 2);                    % By source only
inputtodet = nansum(det,1) + sum(S.dtImp);  % By destination only

% Consumption grid

q0 = bsxfun(@times, S.dc, S.qb' .* S.b');   % consumption by all groups
q0(:, ~islive) = detgroups;                 % "consumption" by detritus
q0(~islive, ~islive) = 0;                   % Fixes detritus ee calc (no NaN)

% Detritus loss to consumption by other groups

deteaten = sum(q0(~islive,:),2)';  

% Calculate BA of detritus (based on detritus fate of detritus groups)
% For now, I'm ignoring this and assuming that all surplus detritus is
% exported from the system, as in my test cases. 

detsurplus = inputtodet - deteaten;         % surplus in each detritus group
ba(~islive) = detsurplus' .* S.df(~islive); % where surplus goes

% Respiration

temp(S.pp < 1)  = 1 - S.pp(S.pp < 1);
temp(S.pp >= 1) = 1;

respiration = S.b .* S.qb - temp .* (S.ee .* S.b .* S.pb + flowtodet(1:S.ngroup));
respiration(S.pp >= 1) = 0;

% EE of detritus

needdetee = ~islive & isnan(S.ee);
tempee(~islive) = (deteaten ./ inputtodet)' - respiration(~islive);
S.ee(needdetee) = tempee(needdetee);

% Export of detritus

detexport = inputtodet' - deteaten' - S.ba(S.nlive+1:end) - respiration(S.nlive+1:end);

%--------------------------
% Other
%--------------------------

migration = S.immig - S.emigRate;
migrationRate = migration ./ S.b;

S.baRate = S.ba./S.b;

%fishMortRate = sum(landing - discard, 2) ./ b;
fishMortRate = bsxfun(@rdivide, S.landing-S.discard, S.b); % fishing rate per biomass by gear

predMort = S.dc * (S.qb .* S.b);                  % total for each group, M2*B in documentation
predMortRate = predMort ./ S.b;                   % rate wrt biomass of prey, M2 in documentation
predMort2 = bsxfun(@times, (S.qb .* S.b)', S.dc); % breakdown for each pred/prey relationship
predMort2 = bsxfun(@rdivide, predMort2, S.b);     % rate wrt biomass of prey

otherMortRate = S.pb .* (1 - S.ee);

q0sum = sum(q0,1); % Note: doesn't include import

searchRate = q0 ./ (S.b * S.b');
searchRate(:,~islive) = 0;

trophic = trophiclevel(S.dc, S.pp, S.nlive, S.ngroup);

%--------------------------
% Print results or assign
% output
%--------------------------

C.ngroup        = S.ngroup;

% Basic estimates

C.trophic       = trophic;
C.areafrac      = S.areafrac;
C.bh            = S.bh;
C.b             = S.b;
C.pb            = S.pb;
C.qb            = S.qb;
C.ee            = S.ee;
C.ge            = S.ge;

% Key indices

C.ba            = S.ba;
C.baRate        = S.baRate;
C.migration     = migration;
C.flowtodet     = flowtodet;

% Mortalities

C.fishMortRate  = fishMortRate;
C.predMortRate  = predMortRate;
C.migrationRate = migrationRate;
C.otherMortRate = otherMortRate;
C.predMort      = predMort2;

% Consumption

C.q0            = q0;
C.q0Sum         = q0sum;

% Respiration

C.respiration   = respiration;

% Search rates

C.searchRate    = searchRate;

% Detritus export

C.detexport     = detexport;
   
if nargout == 0
    displayecopath(S,C);
elseif nargout == 1
    varargout{1} = C;
elseif nargout == 2
    varargout{1} = C;
    varargout{2} = flag;
end

%************************** Subfunctions **********************************

%---------------------------
% Determine whether B and 
% Q/B of all predators of 
% each group are known
%---------------------------

% function isknown = knowpredatorinfo(predIdx, b, qb) 
% isknown = cellfun(@(x) all(~isnan(b(x))) & all(~isnan(qb(x))), predIdx);

%---------------------------
% Calculate total predation 
% loss of each group
%---------------------------

% function m2 = calcpredation(predIdx, predDc, b, qb)
% m2 = nansum(bsxfun(@times, b' .* qb', predDc), 2);
% % m2 = zeros(size(b));
% for igroup = 1:length(predIdx)
% %     m2(igroup) = 
%     
%     m2(igroup) = sum(b(predIdx{igroup}) .* qb(predIdx{igroup}) .* ...
%                      predDc{igroup}');
% end
% % m2 = cellfun(@(x,y) sum(b(x) .* qb(x) .* y), predIdx, predDc);
% m2 = m2';

%---------------------------
% Calculate net migration 
% for each group
%---------------------------

function nm = calcnetmigration(emigRate, emig, immig)

temp = emig > 0 & emigRate == 0;
emigtemp = zeros(size(emig));
emigtemp(temp) = emig(temp);

nm = emigtemp - immig;

%---------------------------
% Determine whether B, P/B, 
% EE of a group's prey 
% groups, and the B and Q/B
% for all predators of the 
% prey group except itself, 
% are known
%---------------------------

% function knowAllPreyInfo = knowpreysotherpredinfo(ngroup, preyIdx, ...
%                                           preysOtherPredIdx, b, qb, pb, ee)
% 
% for igroup = 1:ngroup
%     knowPreysOtherPredInfo{igroup} = ...
%         cellfun(@(x) ~isempty(x) & all(~isnan(b(x))) & all(~isnan(qb(x))), ...
%             preysOtherPredIdx{igroup});
% %         cellfun(@(x) length(x)>0 && ~any(isnan([b(x);qb(x)])), preysOtherPredIdx{igroup});
% 
% end
% knowPreysB = arrayfun(@(x) ~isnan(b(preyIdx{x})), 1:ngroup, 'uni', 0);
% knowPreysPb = arrayfun(@(x) ~isnan(pb(preyIdx{x})), 1:ngroup, 'uni', 0);
% knowPreysEe = arrayfun(@(x) ~isnan(ee(preyIdx{x})), 1:ngroup, 'uni', 0);
% knowAllPreyInfo = cellfun(@(x1,x2,x3,x4) x1'&x2&x3&x4, ...
%                           knowPreysOtherPredInfo, knowPreysB, knowPreysPb, ...
%                           knowPreysEe, 'uni', 0);

%---------------------------
% Change all empty arrays to 
% [] (not 0 x 1 or 1 x 0 
% cell arrays) to avoid
% dimension problems later
%---------------------------

% function c = emptyelement(c)
% isemp = cellfun('isempty', c);
% [c{isemp}] = deal([]);

%---------------------------
% Production-biomass-
% consumption calculations
%---------------------------

function [pb, qb, ge] = pbq(pb, qb, ge, islive)

temp = isnan(pb(islive)) & ~isnan(qb(islive)) & ~isnan(ge(islive));
pb(temp) = ge(temp) .* qb(temp);
        
temp = isnan(qb(islive)) & ~isnan(pb(islive)) & ~isnan(ge(islive));
qb(temp) = pb(temp) ./ ge(temp);
        
temp = ~isnan(qb(islive)) & ~isnan(pb(islive));
ge(temp) = pb(temp) ./ qb(temp);
     
                        
%---------------------------
% Check if values are known
%--------------------------- 

function knowall = checkbasic(b, pb, qb, ee, ge, islive, pp)
knowall = ~any(isnan(b)) && ...
          ~any(isnan(pb(islive))) && ...
          ~any(isnan(qb(pp == 0))) && ...
          ~any(isnan(ee(islive))) && ...
          ~any(isnan(ge(pp == 0)));

                            
                        
                        
         
    

