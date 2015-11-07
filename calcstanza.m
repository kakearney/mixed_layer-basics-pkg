function A = calcstanza(A, varargin)
% CALCSTANZA Calculate B and Q/B values for multi-stanza Ecopath groups
%
% A = calcstanza(In)
% A = calcstanza(In, p1, v1)
%
% This function calculates the B and QB values associated with multistanza
% ecopath groups.
%
% A few caveats...
%
% The values calculated by this script do not exactly replicate the values
% one will get using the EwE 6 "Edit multi-stanza" menu.  For the biomass,
% the difference should be very small, and results I think simply due to
% the fact that my calculations are all carried out in double precision,
% while the EwE code uses a mix of single and integer precision.  Although
% Matlab makes it relatively easy to calculate the direct integrals (over
% mortality rate Z and weight w), for consistency I have instead stuck with
% the convention used in EwE to use a monthly discretization.  In tests,
% the percent error between my results and the EwE ones tend to be <0.01%.
%
% However, the error on the consumption rates is a bit higher (in tests,
% 0.005-0.4%), and I haven't quite figured out the source of this yet.  It
% seems a bit high to be just due to the different precisions, or to the
% fact that I use a much more vectorized approach than the original code.
% While the differences are very small, they may make a difference if one
% plans to generate ensembles using this code and then move back to the EwE
% software.  If you do this, I suggest you double-check the balance of the
% model(s) there.
%
% Input variables:
%
%   In:     Ecopath input structure.  B and QB of non-leading multi-stanza 
%           groups can be unknown (i.e. NaN).  Stanza-related fields
%           (stanza, stanzadata.stanzaID, stanzadata.BABsplit, ageStart,
%           and vbK) must be present.   
%
% Optional input variables (passed as parameter/value pairs)
%
%   plot:   logical scalar, true to plot growth curve details for each
%           stanza group, primarily for debugging purposes.  [false] 
%
%   da:     discretizatin interval (months) for age calculations [1]
%
% Output variables:
%
%   A:      Ecopath input structure, identical to In except that
%           non-leading multi-stanza group B and QB values have been
%           recalculated/filled in.  

% Copyright 2014-2015 Kelly Kearney

% Parse input

ns = length(A.stanzadata.StanzaID);

Opt.plot = false;
Opt.da = 1;

Opt = parsepv(Opt, varargin);

% Set up plotting (if true)

if Opt.plot
    
    h.fig = figure;
    
    vfrac = 0.9;
    hght1 = vfrac./ns;
    hght2 = hght1*0.8;
    
    yax = 0.05:hght1:0.9;
    yax = yax(end:-1:1);
    
    if verLessThan('matlab', 'R2014b')
        h.ax = zeros(ns,1);
    else
        h.ax = gobjects(ns,1);
    end
    for ii = 1:ns
        h.ax(ii) = axes('position', [0.05 yax(ii) 0.9 hght2]);
    end
end
    
% Loop over multi-stanza groups

for is = 1:ns
    
    idx = find(A.stanza == is);
    
    [a, isrt] = sort(A.ageStart(idx));
    idx = idx(isrt);
   
    % Parameters
    
    k = A.vbK(idx(1));               % Curvature parameter
    if iscell(A.stanzadata.BABsplit)
        bab = A.stanzadata.BABsplit{is};
    else
        bab = A.stanzadata.BABsplit(is);
    end
    z = A.pb(idx);
    blead = A.b(idx(end));
    qblead = A.qb(idx(end));
    
    % Fill in non-leading group values
    
    [Out,D] = editstanzacalcs(a, k, bab, blead, qblead, z, Opt.da);
    
    A.b(idx)  = Out.b;
    A.qb(idx) = Out.qb;
    A.ba(idx) = Out.ba;
    
    % Plot to check, mimicking the one in EwE6
    
    if Opt.plot
        
        yfac = max(D.y);
        ynorm = bsxfun(@rdivide, D.y, yfac);
        line(D.x, ynorm, 'parent', h.ax(is));
        line([a a]', (ones(size(a))*[0 1])', 'color', 'k', 'parent', h.ax(is));
        text(max(D.x)*0.98, 0.5, A.stanzadata.StanzaName{is}, ...
            'parent', h.ax(is), 'vert', 'top', 'horiz', 'right', ...
            'fontsize', 8);
        set(h.ax(is), 'xlim', [0 max(D.x)]);     
    end
    
end

% Add labels to plots

if Opt.plot
    xlabel(h.ax(end), 'Age (months)');
    lbl = {...
            'w = von Bertalanffy body weight'
            'l = survivorship'
            'w*l'
            '$\sum_{0}^{a}{Z_a} - a\frac{BA}{B}$'};
    legendflex(h.ax(1), lbl, 'ref', h.ax(1), 'nrow', 1, 'anchor', {'n','s'}, ...
        'xscale', 0.5, 'interpreter', 'Latex', 'buffer', [0 5]);
end
    






