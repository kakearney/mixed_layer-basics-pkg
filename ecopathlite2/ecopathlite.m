function varargout = ecopathlite(S, varargin)
%ECOPATHLITE Rewrite of Ecopath algorithms
%
% ecopathlite(S)
% C = ecopathlite(S)
% [C, flag, fillinfo, sc] = ecopathlite(S)
% [C, CE] = ecopathlite(S, 'x', x, 'idx', idx, ...)
%
% This function reproduces the main calculations performed by the Ecopath
% portion of the EwE model (www.ecopath.org).  
%
% Ecopath is used to calculate a snapshot of an ecosystem, including the
% biomass of all functional groups (living groups, fishing gears, and
% detrital pools) and the fluxes between these groups.  This function is a
% bare-bones version of the algorithm; it is not really meant to be used to
% set up and balance a model for the first time, since it does not provide
% any of the visual checks on the results (e.g. whether EE values are > 1,
% etc); use the original Ecopath software if you're looking for this type
% of behavior.  Nor does it include many of the more complicated setup
% options, such as economic variables (e.g. market prices, fleet dynamics),
% etc.  It simply provides the initial mass-balance calculations, in a
% stripped-down, easy-to-automate format.  
%
% The units below are defined in terms of three variables: M = mass, A =
% area (or volume), and T = time. In the original software, the default is
% M = metric tons wet weight, A = km^2, and T = year.  You can use whatever
% units work best for your own purposes, as long as they remain consistent
% across all variables.  
%
% Input and output variables were chosen to be somewhat consistent with the
% tables seen in Ecopath with Ecosim version 5.  If no output variable is
% provided, the results are printed in the command window, mimicing the
% Basic Estimates table from EwE.
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
% EwE5), and continue to add updates based on the EwE6 source code.
%
% Input variables:
%
%   S:          structure with the following fields.  Values of the fields
%               b, pb, qb, ee, ge, gs, and/or dtImp that are defined as NaN
%               indicate unknown values, which will be filled in by the
%               ecopath algorithm.
% 
%               ngroup:         1 x 1 array, number of functional groups in
%                               the model 
%
%               nlive:          1 x 1 array, number of live (non-detrital)
%                               groups in the model    
%
%               ngear:          1 x 1, number of fishing gear types in the
%                               model 
%
%               areafrac:       ngroup x 1 array, fraction of habitat area
%                               occupied by each group (no units, 0-1) 
%
%               b:              ngroup x 1 array, biomass (M A^-1)
%
%               pb:             ngroup x 1 array, production over biomass
%                               ratios (T^-1)   
%   
%               qb:             ngroup x 1 array, consumption over biomass
%                               ratios (T^-1) 
%
%               ee:             ngroup x 1 array, ecotrophic efficiencies
%                               (no units, 0-1) 
% 
%               ge:             ngroup x 1 array, gross efficiency, i.e.
%                               production over consumption ratio (no
%                               units)
%
%               gs:             ngroup x 1 array, fraction of consumed food
%                               that is not assimilated (no units)  
%
%               dtImp:          ngroup x 1 array, detritus import (should
%                               be zero for all non-detrital groups) (M
%                               A^-1 T^-1)
%
%               bh:             ngroup x 1 array,  habitat biomass, i.e.
%                               biomass per unit habitable area (M A^-1).
%                               This variable is designed as a shortcut if
%                               all your critters are clustered in only a
%                               small portion of the habitat area, such
%                               that bh = b/areafrac.
%
%               pp:             ngroup x 1 array, fraction of diet
%                               consisting of primary production, pp = 2
%                               indicates detritus 
%
%               import:         ngroup x 1 array, fraction of diet coming
%                               from outside the system. 
% 
%               dc:             ngroup x ngroup array, diet composition,
%                               dc(i,j) tells fraction predator j's diet
%                               consisting of prey i (no units, 0-1)
% 
%               df:             ngroup x (ngroup - nlive) array, fraction
%                               of each group that goes to each detrital
%                               group due to other mortality and egestion
%                               (no units, 0-1)
%
%               immig:          ngroup x 1 array, immigration into area  (M
%                               A^-1 T^-1)  
%
%               emig:           ngroup x 1 array, emigration out of area
%                               (M A^-1 T^-1) 
%
%               emigRate:       ngroup x 1 array, emigration per unit
%                               biomass (T^-1)
% 
%               ba:             ngroup x 1 array, biomass accumulation  (M
%                               A^-1 T^-1) 
%
%               baRate          ngroup x 1 array, biomass accumulation per
%                               unit biomass (T^-1) 
%
%               landing         ngroup x ngear array, landings of each
%                               group by each gear type (M A^-1 T^-1) 
%
%               discard         ngroup x ngear array, discards of each
%                               group by each gear type (M A^-1 T^-1) 
%
%               discardFate:    ngear x (ngroup - nlive) array, fraction of
%                               discards from each gear type that go to
%                               each detritus group (no units, 0-1)
%
%               stanzadata:     dataset array holding properties of each
%                               multi-stanza group (only necessary if at
%                               least one multi-stanza group in the model)
%
%               stanza:         ngroup x 1 array, stanza ID corresponding
%                               to each group, 0 if not a multi-stanza
%                               group (only necessary if at least one
%                               multi-stanza group in the model)
%
%               ageStart:       ngroup x 1 array, age (in months) at which
%                               each stanza group starts  (only necessary
%                               if at least one multi-stanza group in the
%                               model)
%
%               vbK:            ngroup x 1 array, K parameter from von
%                               Bertalanffy growth curve.  Should be same
%                               value for all age groups in a single
%                               stanza.  EwE6 uses -1 as a placeholder for
%                               non-stanza groups (only necessary if at
%                               least one multi-stanza group in the model)
%
% Optional input variables (passed as parameter/value pairs):
%
%   x:          1 x 6 cell array of varying parameter values (see
%               createensemble.m).  If included (along with idx input),
%               Ecopath balances will be calculated for each ensemble
%               member, with results returned to the CE output.  Empty by
%               default
%
%   idx:        1 x 6 cell array, indices corresponding to the values in x
%               (see createensemble.m).  Empty by default. 
%
%   debug:      logical scalar, true to track and return extra information
%               for debugging purposes.  Not applicable to multi-ensemble
%               member runs.  False by default.  
%
%   skipextra:  logical scalar. If true, the script will fill in missing
%               values but not perform any additional calculations, and
%               will return only the b, pb, qb, ee, and ge fields of C. 
%
%   silent:     logical scalar.  If true, all warning messages are
%               suppressed.
%
% Output variables:
%
%   C:          structure with the following fields:
%
%               trophic:        ngroup x 1 array, trophic level of each
%                               group (no unit)  
%
%               areafrac:       ngroup x 1 array, fraction of total area
%                               occupied by group (no units, 0-1) 
%
%               bh:             ngroup x 1 array, biomass in habitable area
%                               (M A^-1) 
%
%               b:              ngroup x 1 array, total biomass (M A^-1)
%
%               pb:             ngroup x 1 array, production/biomass ratio
%                               (T^-1) 
%   
%               qb:             ngroup x 1 array, consumption/biomass ratio
%                               (T^-1) 
%
%               ee:             ngroup x 1 array, ecotrophic efficiency (no
%                               unit) 
%
%               ge:             ngroup x 1 array, growth efficiency, i.e.
%                               production/consumption (no unit) 
%
%               ba:             ngroup x 1 array, biomass accumulation  (M
%                               A^-1 T^-1) 
%
%               baRate:         ngroup x 1 array, biomass accumulation per
%                               unit biomass (T^-1)  
%
%               migration:      ngroup x 1 array, net migration (M A^-1
%                               T^-1) 
%
%               flowtodet:      (ngroup + ngear) x 1 array, flow to
%                               detritus from each group and each gear type
%                               (M A^-1 T^-1)  
%
%               fishMortRate:   ngroup x 1 array, mortality per unit
%                               biomass due to fishing (T^-1) 
%
%               predMortRate:   ngroup x 1 array, mortality per unit
%                               biomass due to predation, M2 in some
%                               documentation (T^-1)   
%
%               migrationRate:  ngroup x 1 array, net migration per unit
%                               biomass (T^-1) 
%
%               otherMortRate:  ngroup x 1 array, mortality per unit
%                               biomass due to anything else (T^-1) 
%
%               predMort:       ngroup x ngroup array, predation mortality
%                               broken down by predator and prey,
%                               C.predMoreRate = sum(C.predMort,2). (T^-1)
%
%               q0:             ngroup x ngroup array, q0(i,j) is the flux
%                               of biomass from group i to group j.  For j
%                               = live, this is due to predation.  For j =
%                               detrital, this includes non-predatory
%                               mortality and egestion.  (M A^-1 T^-1)
%
%               q0Sum:          ngroup x 1 array, total consumption by each
%                               predator (M A^-1 T^-1) 
%
%               respiration:    ngroup x 1 array, respiration (M A^-1 T^-1)
%
%               searchRate:     ngroup x ngroup array, search rates of each
%                               predator for each prey, assuming simple
%                               linear dynamics, i.e. Qij = a*Bi*Bj (A M^-1
%                               T^-1)  
%
%               detexport:      (ngroup-nlive) x 1 array, amount of
%                               detritus exported from the system (M A^-1
%                               T^-1) 
%
%               omnivory:       ngroup x 1, omnivory index, i.e. variance
%                               of prey's trophic level (Note: not quite
%                               right for partial primary producers right
%                               now) (no unit)
%
%               neteff:         ngroup x 1 array, net food conversion
%                               efficiency (no unit) 
%
%               import:         (ngroup + ngear) x 1 array, import into the
%                               system to each group/gear.  Includes
%                               consumption of prey from outside and
%                               detritus import (M A^-1 T^-1) 
%
%               export:         (ngroup + ngear) x 1 array, export out of
%                               the system by each group/gear.  Includes
%                               nonpredatory and egested loss not directed
%                               to detritus, fisheries landings, and
%                               detritus export. 
%
%               flow:           (ngroup + ngear + 1) x (ngroup + ngear + 2)
%                               array, matrix of all flows in the system,
%                               with  rows corresponding to source groups
%                               and columns to sink groups.  Rows/columns
%                               are in the order of live groups, detrital
%                               groups, fishing gears, outside the system,
%                               and respiration-destination.
%
%               Idx:            structure of indices corresponding to the
%                               rows and columns of the flow array,
%                               categorized as living groups (liv),
%                               detrital groups(det), fishing gears (gear),
%                               outside the system (out), and the
%                               mysterious beyond where respiration goes
%                               (res).                       
%
%   CE:         1 x nens structure with same format as C. Multi-ensemble
%               runs only, where nens is the number of ensemble members
%               described by the parameters in input x.
%
%   flag:       false if all unknowns are filled, true if not.  This output
%               was primarily added for testing purposes.  Unfilled
%               unknowns are usually due to incorrect input, but may point
%               to a bug in this implementation of the Ecopath algorithm.
%               Output not available with multi-ensemble calculation
%               (Assuming the main model of the ensemble can be filled, the
%               other ensemble members should too, with one exception: the
%               parameter-variation process can push groups with
%               cannibalism over the threshold where predation mortality
%               exceeds cannibalism, leading to some ensemble members with
%               missing values... unless running in silent mode, these will
%               be noted by a warning message printed to the screen).
%
%   fillinfo:   dataset array indicating which algorithm (see Appendix 4 
%               of the EwE User's Manual) is used to fill in each value,
%               and on which iteration it was filled.  Also primarily for
%               testing purposes. Output not available with multi-ensemble
%               calculation. 
%
%   sc:         Sanity check calculation that lists the main terms of the
%               Ecopath equation for each group.  Columns correspond to
%               Bi*PBi*EEi, B1*QB1*DC1i, B2*QB2*DCi2, ... Yi, Ei, BAi.  The
%               last column is the sum of each row, and should sum to 0 (or
%               very close, near machine precision).  Output not available
%               with multi-ensemble calculation. 

% Copyright 2012-2015 Kelly Kearney
% kakearney@gmail.com

Opt.x = [];
Opt.idx = [];
Opt.skipextra = false;
Opt.debug = false;
Opt.silent = false;

Opt = parsepv(Opt, varargin);

multiflag = ~isempty(Opt.x);
if ~multiflag && nargout > 1
    Opt.debug = true;
end

%------------------------------
% Setup
%------------------------------

if ~isscalar(S)
    error('Input structure must be scalar');
end

if Opt.silent
    S = ecopathinputcheck(S, true);
else
    S = ecopathinputcheck(S);
end

if multiflag
    Ens = subecopathens(S, Opt.x, Opt.idx);
    nens = length(Ens);
end

Sorig = S;

%------------------------------
% Setup calculations
%------------------------------

islive = (1:S.ngroup)' <= S.nlive; % Logical mask for live groups

% If emigration and biomass is given as an input rather than emigration per
% unit biomass (i.e. emigration rate), calculate the rate, and vice versa.
% Same for BA.

S = convertrates(S, islive);

% Calculate growth efficiency (i.e. production/consumption ratio),
% production/biomass ratio, and consuption/biomass ratio if needed and
% possible

S = pbq(S, islive);

% Calculate total catches for each group

S.catches = sum(S.landing + S.discard, 2);

if multiflag
    
    for ie = 1:nens
        Ens(ie) = convertrates(Ens(ie), islive);
        Ens(ie) = pbq(Ens(ie), islive);
        Ens(ie).catches = sum(Ens(ie).landing + Ens(ie).discard, 2);
    end
end

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

if Opt.debug
    count = 0;
    filliter = nan(S.ngroup, 5);
    fillalgo = nan(S.ngroup, 5);
    status = [S.b S.pb S.qb S.ee S.ge];
end

Fail.flag = 0;
Fail.idx = 0;
Fail = struct('flag', 0, 'idx', 0);
if multiflag
    FailE = struct('flag', num2cell(zeros(nens,1)), 'idx', cell(nens,1));
end

while ~checkbasic(S.b, S.pb, S.qb, S.ee, S.ge, islive, S.pp)

    if Opt.debug
        count = count + 1;
    end
    
    param = [S.b S.pb S.qb S.ee S.ge];
    
    % Run the p-b-q algebra again (in case new values have been filled in),
    % and the rate-to-total calcs (which I'm pretty sure won't change,
    % since I don't think EwE6 allows you to enter a rate if the B is
    % missing, but just in case...)
    
    S = pbq(S, islive);
    S = convertrates(S, islive);
    
    if multiflag
        for ie = 1:nens
            Ens(ie) = pbq(Ens(ie), islive);
            Ens(ie) = convertrates(Ens(ie), islive);
        end
    end

    % Total export
    
    S.ex = S.catches + S.emig - S.immig + S.ba; 
    if multiflag
        for ie = 1:nens
            Ens(ie).ex = Ens(ie).catches + Ens(ie).emig - Ens(ie).immig + Ens(ie).ba; 
        end
    end
    
    %-----------------------
    % Algorithm 1: 
    % Estimation of P/B
    %-----------------------

    knowPredInfo = ~any(bsxfun(@and, ispred, isnan(S.b))) & ...
                   ~any(bsxfun(@and, ispred, isnan(S.qb)));
               
    canRun = (islive & ...          % only live groups           
              isnan(S.pb) & ...     % P/B unknown
              ~isnan(S.b) & ...     % B known
              ~isnan(S.ee) & ...    % EE known
              knowPredInfo');       % know B and Q/B for all group's predators
          
    S = alg1(S, canRun);
    if multiflag
        for ie = 1:nens
            Ens(ie) = alg1(Ens(ie), canRun);
        end
    end      
    
    if Opt.debug
        ischanged = ~arrayfun(@(x,y) isequaln(x,y), status, [S.b S.pb S.qb S.ee S.ge]);
        status = [S.b S.pb S.qb S.ee S.ge];
        fillalgo(ischanged) = 1;
        filliter(ischanged) = count;
    end
    
    %-----------------------
    % Algorithm 2: 
    % Estimation of EE
    %-----------------------
    
    knowPredInfo = ~any(bsxfun(@and, ispred, isnan(S.b))) & ...
                   ~any(bsxfun(@and, ispred, isnan(S.qb)));
    
    canRun = (islive & ...          % only live groups           
              isnan(S.ee) & ...     % EE unknown
              ~isnan(S.b) & ...     % B known
              ~isnan(S.pb) & ...    % P/B known
              knowPredInfo');       % know B and Q/B for all group's predators
          
    S = alg2(S, canRun);
    if multiflag
        for ie = 1:nens
            Ens(ie) = alg2(Ens(ie), canRun);
        end
    end  
    
    if Opt.debug
        ischanged = ~arrayfun(@(x,y) isequaln(x,y), status, [S.b S.pb S.qb S.ee S.ge]);
        status = [S.b S.pb S.qb S.ee S.ge];
        fillalgo(ischanged) = 2;
        filliter(ischanged) = count;
    end
    
    %-----------------------
    % Algorithm 3: Dealing 
    % with B and Q/B as 
    % unknowns  
    %-----------------------

    % It's never stated in the docs, but the prey group k for which B, PB,
    % EE, and predator info needs to be know cannot be a detrital group (or
    % at least, that seems to lead to incorrect results). To date, I've
    % never found model that needs this algorithm, so it remains untested.

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
              any(knowAllPreyInfo(islive,:), 1)');  % know B, P/B, EE, and pred info (except group) for at least one live prey    

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

    S = alg3(S, group1, group2, k, idx);
    if multiflag
        for ie = 1:nens
            Ens(ie) = alg3(Ens(ie), group1, group2, k, idx);
        end
    end
               
    if Opt.debug
        ischanged = ~arrayfun(@(x,y) isequaln(x,y), status, [S.b S.pb S.qb S.ee S.ge]);
        status = [S.b S.pb S.qb S.ee S.ge];
        fillalgo(ischanged) = 3;
        filliter(ischanged) = count;
    end

    %-----------------------
    % Algorithm 4: 
    % Estimating biomasses 
    % only
    %-----------------------

    knowPredInfo = ~any(bsxfun(@and, isPredNoCannib, isnan(S.b))) & ...
                   ~any(bsxfun(@and, isPredNoCannib, isnan(S.qb)));


    canRun = (islive & ...          % only live groups
              isnan(S.b) & ...      % B unknown
              ~isnan(S.pb) & ...    % P/B known
              ~isnan(S.ee) & ...    % EE known
              ~isnan(S.qb) & ...    % Q/B known
              knowPredInfo');       % B and Q/B of predators except itself known

     if ~Fail.flag     
        [S, Fail] = alg4(S, canRun);
     end
    if multiflag
        for ie = 1:nens
            if ~FailE(ie).flag
                [Ens(ie), FailE(ie)] = alg4(Ens(ie), canRun);
            end
        end
    end
    
    if Opt.debug
        ischanged = ~arrayfun(@(x,y) isequaln(x,y), status, [S.b S.pb S.qb S.ee S.ge]);
        status = [S.b S.pb S.qb S.ee S.ge];
        fillalgo(ischanged) = 4;
        filliter(ischanged) = count;
    end
    
    % Only try to fill B and QB if we've done all we can with the first 4
    % algorithms
    
    nochange = isequaln(param, [S.b S.pb S.qb S.ee S.ge]);
    if nochange
        
        %-----------------------
        % Algorithm 5: The 
        % generalized inverse.
        %-----------------------
        
        if ~Fail.flag
            [S, Fail] = alg5(S, islive);
        end
        
        if multiflag
            for ie = 1:nens
                if ~FailE(ie).flag
                    [Ens(ie), FailE(ie)] = alg5(Ens(ie), islive);
                end
            end
        end
        
        % Check again

        if Opt.debug
            ischanged = ~arrayfun(@(x,y) isequaln(x,y), status, [S.b S.pb S.qb S.ee S.ge]);
            status = [S.b S.pb S.qb S.ee S.ge];
            fillalgo(ischanged) = 5;
            filliter(ischanged) = count;
        end
        
    end
    
    % If we made it through all algorithms without filling in any new
    % numbers, break out of loop
    
    nochange = isequaln(param, [S.b S.pb S.qb S.ee S.ge]);
    
    if nochange
        Fail.flag = 4;
    end
    
    if (~multiflag && Fail.flag) ||  (multiflag && Fail.flag && all([FailE.flag]))
        break
    end
    
end

% Issue warnings if any values failed to fill

if ~Opt.silent
    if multiflag
        if any([Fail.flag FailE.flag])
            ecopathwarn([Fail; FailE]);
        end
    else
        ecopathwarn(Fail);
    end
end

% Detritus calculations

S = detcalcs(S, islive);
if multiflag
    
    newfields = setdiff(fieldnames(S), fieldnames(Ens));
    for in = 1:length(newfields)
        [Ens.(newfields{in})] = deal([]);
    end
    for ie = 1:nens
        Ens(ie) = detcalcs(Ens(ie), islive);
    end
end

%-----------------------
% Calculate various
% Ecopath outputs
%-----------------------

C = basicestimates(S, Opt.skipextra);
if multiflag
    for ie = nens:-1:1
        CE(ie) = basicestimates(Ens(ie), Opt.skipextra);
    end
end

if ~Opt.skipextra
    C = keyindices(S,C);
    C = mortalities(S,C);
    C = otheroutput(S,C,islive);
    
    if multiflag
        newfields = setdiff(fieldnames(C), fieldnames(CE));
        for in = 1:length(newfields)
            [CE.(newfields{in})] = deal([]);
        end
        CE = orderfields(CE, C);
        for ie = 1:nens
            CE(ie) = keyindices(Ens(ie),CE(ie));
            CE(ie) = mortalities(Ens(ie),CE(ie));
            CE(ie) = otheroutput(Ens(ie),CE(ie),islive);
        end
    end
end
  

% Assign outputs

if multiflag
    tmp = {C, CE};
    if nargout > 2
        warning('Multi-ensemble mode does not support debugging outputs');
        varargout(1:2) = tmp;
        [varargout{3:nargout}] = deal(NaN);
    else
        varargout = tmp(1:nargout);
    end
else

    if nargout == 0
        displayecopath(Sorig,C);
    elseif nargout == 1
        varargout{1} = C;
    elseif nargout > 1

        isfilled = ~isnan(fillalgo);
        [ig,ivar] = find(isfilled);
        var = {'B','PB','QB','EE','GE'}';
        fillinfo = dataset({var(ivar), 'Variable' }, ...
                           {ig, 'Group'}, ...
                           {fillalgo(isfilled), 'Algorithm'}, ...
                           {filliter(isfilled), 'Iteration'});

        tmp = {C, Fail.flag, fillinfo, sanitycheck(S)};
        varargout = tmp(1:nargout);
    end
end

%************************** Subfunctions **********************************

%---------------------------
% Production-biomass-
% consumption calculations
%---------------------------

function S = pbq(S, islive)

temp = isnan(S.pb(islive)) & ~isnan(S.qb(islive)) & ~isnan(S.ge(islive));
S.pb(temp) = S.ge(temp) .* S.qb(temp);
        
temp = isnan(S.qb(islive)) & ~isnan(S.pb(islive)) & ~isnan(S.ge(islive));
S.qb(temp) = S.pb(temp) ./ S.ge(temp);
        
temp = ~isnan(S.qb(islive)) & ~isnan(S.pb(islive));
S.ge(temp) = S.pb(temp) ./ S.qb(temp);
     
                        
%---------------------------
% Check if values are known
%--------------------------- 

function knowall = checkbasic(b, pb, qb, ee, ge, islive, pp)
knowall = ~any(isnan(b)) && ...
          ~any(isnan(pb(islive))) && ...
          ~any(isnan(qb(pp == 0))) && ...
          ~any(isnan(ee(islive))) && ...
          ~any(isnan(ge(pp == 0)));

%---------------------------
% Rate-to-total conversions
%---------------------------  
      
function S = convertrates(S, islive)

% NOTE: In previous versions of this code, I seemed to be using reverse
% terminology for emig vs emigRate (calcs were right, though)... not sure
% whether that was a mistake on my part or a convention I stole from the
% EwE6 code, but either way it was really confusing me, so I've switched
% back.

if any(isnan(S.b) & (S.emigRate > 0 | S.baRate > 0))
    warning('Missing b combined with assigned Emigration and/or BA rate: This scenario may not work... check');
end

e2er   = islive & (S.emig > 0) & ~isnan(S.b) & (S.emigRate == 0);
er2e   = islive & ~isnan(S.b) & (S.emigRate > 0) & (S.emig == 0);
ba2bar = islive & (S.ba ~= 0) & ~isnan(S.b) & (S.baRate == 0);
bar2ba = islive & ~isnan(S.b) & (S.baRate ~= 0) & (S.ba == 0);

S.emigRate(e2er) = S.emig(e2er) ./ S.b(e2er);
S.emig(er2e) = S.emigRate(er2e) .* S.b(er2e);
S.baRate(ba2bar) = S.ba(ba2bar) ./ S.b(ba2bar);
S.ba(bar2ba) = S.baRate(bar2ba) .* S.b(bar2ba);
                            
%---------------------------
% Sanity check: are things
% balancing properly?
%---------------------------           

function allterms = sanitycheck(S)

qtmp = bsxfun(@times, S.b' .* S.qb', S.dc);
qtmp(S.dc == 0) = 0; 
catches = sum(S.landing + S.discard, 2);

% The master Ecopath equation
% Bi*PBi*EEi - sum_over_j(Bj*QBj*DCij) - Yi - Ei - BAi = 0
%
% Last displayed column should be all 0s

allterms = [S.b.*S.pb.*S.ee -qtmp -catches -S.emig S.immig -S.ba];
allterms = [allterms sum(allterms,2)];

%---------------------------
% Ecopath algebra algorithms
%---------------------------  

% 1: Estimation of P/B

function S = alg1(S, canRun)
m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2);   
S.pb(canRun) = (S.ex(canRun) + m2(canRun))./(S.b(canRun) .* S.ee(canRun)); 

% 2: Estimation of EE

function S = alg2(S, canRun)
m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2); 
S.ee(canRun) = (S.ex(canRun) + m2(canRun))./(S.b(canRun) .* S.pb(canRun));

% 3: B and Q/B as unknowns

function S = alg3(S, group1, group2, k, idx)
partm2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc .* ~eye(S.ngroup)), 2); % predation, not including cannibalism      
m2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc), 2); % all predation

dcii = diag(S.dc);
dcik = S.dc(idx);

S.b(group1) = partm2(group1) + S.ex(group1) + ...
              dcii(group1) .*(S.b(k(group1)) .* S.pb(k(group1)) .* ...
              S.b(k(group1)) .* S.ee(k(group1)) - S.ex(k(group1)) - ...
              m2(k(group1)))./dcik(group1);

S.qb(group2) = (S.b(k(group2)) .* S.pb(k(group2)) .* S.b(k(group2)) .* ...
               S.ee(k(group2)) - S.ex(k(group2)) - m2(k(group2))) ./ ...
               (dcik(group2) ./ S.b(group2));

% 4: Biomasses only

function [S, Fail] = alg4(S, canRun)  

partm2 = nansum(bsxfun(@times, S.b' .* S.qb', S.dc .* ~eye(S.ngroup)), 2); % predation, not including cannibalism      
dcCannib = diag(S.dc);

cannibCheck1 = canRun & (S.pb .* S.ee) == (S.qb .* dcCannib);
cannibCheck2 = canRun & (S.pb .* S.ee) < (S.qb .* dcCannib);

if any(cannibCheck1)
    Fail.flag = 1;
    Fail.idx = find(cannibCheck1);
elseif any(cannibCheck2)
    Fail.flag = 2;
    Fail.idx = find(cannibCheck2);
else
    Fail.flag = 0;
    Fail.idx = [];
end

if ~Fail.flag
    S.b(canRun) = (S.ex(canRun) + partm2(canRun)) ./ ...
                  (S.pb(canRun) .* S.ee(canRun) - S.qb(canRun) .* dcCannib(canRun));
end

% 5: The generalized inverse

function [S, Fail] = alg5(S, islive)

bmiss = isnan(S.b);
qbmiss = isnan(S.qb);

if any(bmiss & qbmiss)
    Fail.flag = 3;
    Fail.idx = find(bmiss & qbmiss);
else
    Fail.flag = 0;
    Fail.idx = [];
end

if ~Fail.flag

    % Set up matrices: AX = Q, from
    % Bi*(P/B)i*EEi - sum_over_j(Bj*(Q/B)j*DCij) - Yi - Ei - BAi = 0

    rhs = [bsxfun(@times, S.b' .* S.qb', S.dc), ...
           -S.b.*S.pb.*S.ee, ...
           S.ex];

    rhs(:,bmiss | qbmiss) = 0; 
    rhs(bmiss,S.ngroup+1) = 0; 

    Q = sum(rhs,2);

    lhs1 = S.pb .* S.ee;
    lhs1(~bmiss) = 0;
    lhs1 = diag(lhs1);

    lhs2 = -bsxfun(@times, S.qb', S.dc);
    lhs2(:,~bmiss) = 0;

    lhs3 =  -bsxfun(@times, S.b', S.dc);
    lhs3(:,~qbmiss) = 0;

    A = lhs1 + lhs2 + lhs3;

    % Drop unecessary columns and rows w/ NaNs (from missing PB or EE)

    isn = isnan(S.pb) | isnan(S.ee);
    hasterm = any(A,2);

    keeprow = ~isn & islive;% & hasterm; 

    A = A(keeprow, bmiss | qbmiss);
    Q = Q(keeprow);

    % Solve

    X = A\Q;

    % Resubstitute X to B or QB

    borqb = nan(S.ngroup);
    borqb(bmiss | qbmiss) = X;

    S.b(bmiss) = borqb(bmiss);
    S.qb(qbmiss) = borqb(qbmiss);
end

%--------------------------
% Detritus calculations
%--------------------------

function S = detcalcs(S, islive)

% Fill in biomass-in-habitat area if necessary, and correct some values

needBh = isnan(S.bh) & ~isnan(S.b);
S.bh(needBh) = S.b(needBh) ./ S.areafrac(needBh);

S.pb(S.pp == 2) = 0;    % No production for detritus
S.qb(S.pp >= 1) = 0;    % No consumption for primary producers
S.ge(S.pp >= 1) = 0;    % Q = 0 for primary producers, so P/Q = Inf, set to 0 instead just as placeholder
S.gs(S.pp >= 1) = 0;    % Q = 0 so unassim irrelevant for detritus and producers

% Matrix of flows: 
% rows = sources: groups, gears, outside
% cols = sinks:   groups, gears, outside, resp

q0 = bsxfun(@times, S.dc, S.qb' .* S.b'); % Consumption of prey

mort  = bsxfun(@times, S.b .* S.pb .* (1 - S.ee), [S.df 1-sum(S.df,2)]); 
egest = bsxfun(@times, S.b .* S.qb .* S.gs, [S.df 1-sum(S.df,2)]); 
discards = bsxfun(@times, sum(S.discard,1)', [S.discardFate 1-sum(S.discardFate,2)]);
imports = S.qb .* S.b .* S.import;

respiration = S.qb.*S.b - S.pb.*S.b - (S.gs.*S.qb.*S.b);
respiration(S.pp >= 1) = 0;

Idx.liv = 1:S.nlive;
Idx.det = (S.nlive+1):S.ngroup;
Idx.gear = (1:S.ngear) + S.ngroup;
Idx.out = S.ngroup+S.ngear+1; % Import/export
Idx.res = S.ngroup+S.ngear+2; % Respiration
Idx.lan = S.ngroup+S.ngear+3; % Fisheries landings
Idx.gpp = S.ngroup+S.ngear+4; % Primary production

flows = zeros(S.ngroup+S.ngear+4,S.ngroup+S.ngear+4);   
flows([Idx.liv Idx.det], [Idx.liv Idx.det]) = q0;                    % Pred eats prey
flows([Idx.liv Idx.det], [Idx.det Idx.out]) = mort + egest;          % Non-pred loss and egestion goes to detritus or export
flows([Idx.liv Idx.det],  Idx.gear        ) = S.landing + S.discard; % Catches and discards go from groups to gears
flows( Idx.gear,         [Idx.det Idx.out]) = discards;              % Discards are then sent to detritus
flows( Idx.gear,          Idx.lan         ) = sum(S.landing,1);      % Landings are exported from the system 
flows([Idx.liv Idx.det],  Idx.res         ) = respiration;           % Remaining production goes to respiration
flows( Idx.out,          [Idx.liv Idx.det]) = imports;               % Some groups feed outside system
flows( Idx.out,           Idx.det         ) = S.dtImp(Idx.det);      % And detrital groups can also import

primprod = (S.pb.*S.b.*S.pp)';
flows( Idx.gpp, Idx.liv)                    = primprod(Idx.liv);                        % Primary production comes from outside
% TODO: source minus sink is unbalanced for partial primary producers...
% did I miss something? A diet adjustment perhaps?  I've tried a few
% different things to try to eliminate this imbalance, but those lead to a
% mismatch with EwE6, which seems to ignore partial primary production when
% calculating consumption rates.

% What's the remaining balance in the detrital groups?

surplus = nansum(flows(:,Idx.det),1)' - nansum(flows(Idx.det,:),2);

% Surplus either goes to detrital groups or to export

surplusfate = bsxfun(@times, surplus, [S.df(~islive,:) 1-sum(S.df(~islive,:),2)]);

flows(Idx.det, [Idx.det Idx.out]) = surplusfate;

% Move detrital flow to self to BA

toself = diag(flows);  
detba = toself(Idx.det);
S.ba(Idx.det) = detba;  % Haven't thoroughly tested this one...

% EE of detritus

detin = flows(:,Idx.det);
deteaten = flows(Idx.det,Idx.liv);
detresp = flows(Idx.det, Idx.res);
detee = sum(deteaten,2)'./(sum(detin,1) - sum(detresp,2)');

S.ee(~islive) = detee;

% Add the flows field to S for use later

S.flows = flows;
S.Idx = Idx;

% 
% 
% 
% 
% % Detritus produced from mortality and egestion
% 
% mort = bsxfun(@times, S.b .* S.pb .* (1 - S.ee), [S.df 1-sum(S.df,2)]);   
% egest = bsxfun(@times, S.b .* S.qb .* S.gs, [S.df 1-sum(S.df,2)]);   
% detgroups = mort + egest;
% 
% % Detritus produced from fisheries discards
% 
% detfisheries = bsxfun(@times, sum(S.discard,1)', [S.discardFate 1-sum(S.discardFate,2)]);
% 
% % Separate detritus flows from export
% 
% exports = [detgroups(:,end); detfisheries(:,end)];
% detgroups = detgroups(:,1:end-1);
% detfisheries = detfisheries(:,1:end-1);
% 
% % Input to detritus groups from groups, fleets, and import
% 
% det = [detgroups; detfisheries];            % By source and destination
% flowtodet = sum(det, 2);                    % By source only
% inputtodet = nansum(det,1) + sum(S.dtImp);  % By destination only
% 
% % Consumption grid
% 
% q0 = bsxfun(@times, S.dc, S.qb' .* S.b');   % consumption by all groups
% q0(:, ~islive) = detgroups;                 % "consumption" by detritus
% q0(~islive, ~islive) = 0;                   % Fixes detritus ee calc (no NaN)
% 
% % Detritus loss to consumption by other groups
% 
% deteaten = sum(q0(~islive,:),2)';  
% 
% % Respiration
% 
% % temp = zeros(size(S.pp));
% % temp(S.pp < 1)  = 1 - S.pp(S.pp < 1);
% % temp(S.pp >= 1) = 1;
% % 
% % respiration = S.b .* S.qb - temp .* (S.ee .* S.b .* S.pb + flowtodet(1:S.ngroup));
% % respiration(S.pp >= 1) = 0;
% 
% % Respiration (v6.3.1 changed the formula)
% 
% respiration = S.qb.*S.b - S.pb.*S.b - (S.gs.*S.qb.*S.b);
% respiration(S.pp >= 1) = 0;
% 
% % Fate of detritus: surplus detritus (i.e. not eaten) goes either to other
% % detritus groups, to self (as biomass accumulation), or is exported from
% % the system.
% 
% % Note: Inclusion of respiration confuses me... isn't that always 0 for
% % detritus?  Also, in CalcBAofDetritus, they redefine Surplus (my
% % detsurplus) as inputtodet - deteaten - fCatch, to account for a model
% % that included catch of a detritus group, later discarded to a different
% % detritus group.  Not sure why they didn't make that change in
% % CalcFateOfDetritus too; for now I'm leaving it out (seems like an odd
% % edge case anyway).
% 
% detsurplus = inputtodet - deteaten - respiration(~islive)';   % surplus in each detritus group
% surplusfate = bsxfun(@times, S.df(~islive,:), detsurplus');
% 
% isself = eye(size(surplusfate));
% surplusfateself = diag(surplusfate);
% surplusfate     = surplusfate .* ~isself; % set diagnonal to 0
% 
% inputtodet = inputtodet + sum(surplusfate, 1);
% detpassedon = sum(surplusfate,2);
% 
% det(~islive,:) = surplusfate;
% flowtodet = sum(det, 2); % flowtodet(~islive) = detpassedon;
% q0(~islive,~islive) = surplusfate;
% 
% % surplus = inputtodet - deteaten - fcatch(~islive)'; % Mostly same as detsurplus above, but apparently some models inlcude "catch" of detritus, redirected to other detritus groups
% % S.ba(~islive) = surplus' .* S.df(~islive); % where surplus goes
% S.ba(~islive) = surplusfateself;
% 
% % EE of detritus
% 
% needdetee = ~islive & isnan(S.ee);
% tempee(~islive) = (deteaten ./ inputtodet)' - respiration(~islive);
% S.ee(needdetee) = tempee(needdetee);
% 
% % Export of detritus
% 
% hasexport = sum(S.df(~islive,:),2) < 1;
% baDet = S.ba(~islive);
% respDet = respiration(~islive);
% detexport = zeros(S.ngroup-S.nlive,1);
% detexport(hasexport) = inputtodet(hasexport)' - deteaten(hasexport)' - ...
%     baDet(hasexport) - detpassedon(hasexport) - respDet(hasexport);
% 
% exports(~islive) = detexport;
% 
% % Add new variables to S
% 
% S.flowtodet = flowtodet;
% S.detexport = detexport;
% S.respiration = respiration;
% S.q0 = q0;
% S.export = exports;


%-----------------------
% Basic estimates
%-----------------------

function C = basicestimates(S, flag)

if flag
    trophic = nan(S.ngroup,1);
else
    trophic = trophiclevel(S.dc, S.pp, S.nlive, S.ngroup);
end

C.ngroup        = S.ngroup;
C.trophic       = trophic;
C.areafrac      = S.areafrac;
C.bh            = S.bh;
C.b             = S.b;
C.pb            = S.pb;
C.qb            = S.qb;
C.ee            = S.ee;
C.ge            = S.ge;

%-----------------------
% Key indices
%-----------------------

function C = keyindices(S,C)

C.ba            = S.ba;
C.baRate        = S.ba./S.b;
C.migration     = S.emig - S.immig;
C.flowtodet     = sum(S.flows([S.Idx.liv S.Idx.det S.Idx.gear], S.Idx.det), 2);

% C.flowtodet     = S.flowtodet;

bqb = sum(bsxfun(@times, C.trophic, S.dc), 1);
C.omnivory = sum(bsxfun(@minus, C.trophic, bqb).^2 .* S.dc, 1)';

C.neteff = C.pb./(C.qb .* (1-S.gs));

%-----------------------
% Mortalities
%-----------------------

function C = mortalities(S,C)

C.fishMortRate = bsxfun(@rdivide, S.catches, S.b); % fishing rate per biomass by gear

predMort = S.dc * (S.qb .* S.b);                  % total for each group, M2*B in documentation
predMortRate = predMort ./ S.b;                   % rate wrt biomass of prey, M2 in documentation
predMort2 = bsxfun(@times, (S.qb .* S.b)', S.dc); % breakdown for each pred/prey relationship
predMort2 = bsxfun(@rdivide, predMort2, S.b);     % rate wrt biomass of prey

C.predMortRate  = predMortRate;
C.predMort      = predMort2;

C.migrationRate = C.migration ./ S.b;
C.otherMortRate = S.pb .* (1 - S.ee);

%-----------------------
% Consumption, search
% rate, respiration
%-----------------------

function C = otheroutput(S,C,islive)

C.q0 = S.flows([S.Idx.liv S.Idx.det], [S.Idx.liv S.Idx.det]);
C.q0sum = sum(C.q0,1);

C.respiration = S.flows([S.Idx.liv S.Idx.det], S.Idx.res);

C.detexport = S.flows(S.Idx.det, S.Idx.out);

C.searchRate = C.q0 ./ (S.b * S.b');
C.searchRate(:,~islive) = 0;

C.import = S.flows(S.Idx.out, [S.Idx.liv S.Idx.det S.Idx.gear])';
C.export = S.flows([S.Idx.liv S.Idx.det S.Idx.gear], S.Idx.out);

C.flow = S.flows;
C.Idx = S.Idx;

% S.detexport = detexport;
% S.respiration = respiration;
% S.q0 = q0;
% S.export = exports;

% C.q0 = S.q0;
% C.q0Sum = sum(S.q0,1);
% 
% C.respiration = S.respiration;
% 
% C.detexport = S.detexport;
% 
% if isfield(S, 'import')
%     C.impconsump = S.qb .* S.b .* S.import;
% else
%     C.impconsump = zeros(size(S.qb));
% end
% C.export = S.export;

%-----------------------
% Warnings
%-----------------------

function ecopathwarn(F)

for ii = 1:length(F)
    if ii == 1
        ensstr = 'center model';
    else
        ensstr = sprintf('ensemble #%d', ii-1);
    end
    
    switch F(ii).flag
        case 1
            idxStr = sprintf('%d,', F(ii).idx);
            idxStr = idxStr(1:end-1);
            
            warning('EWE:cannibalWithoutB', ...
                'Exiting without solution for %s\nGroup(s) (%s) are missing biomass but are only preyed on by themselves; group(s) must be split in two to solve', ensstr, idxStr);
        case 2
            idxStr = sprintf('%d,', F(ii).idx);
            idxStr = idxStr(1:end-1);
            
            warning('EWE:cannibalTooHigh', ...
                'Exiting without solution for %s\nGroup(s) (%s) have cannibalism losses that exceed predation mortality', ensstr, idxStr);             
        case 3
            idxStr = sprintf('%d,', F(ii).idx);
            idxStr = idxStr(1:end-1);
            
            warning('EWE:B_QB_missing', ...
                'Exiting without solution for %s\nMissing B and QB for group(s) (%s); cannot solve', ensstr, idxStr);
        case 4
            warning('EWE:unknownsremain', ...
                'Exiting without solution for %s\nUnable to fill in all unknowns; check input', ensstr);
        
    end
end
    
