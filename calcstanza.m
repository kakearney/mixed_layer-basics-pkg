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
    h = plotgrid('setup', cell(ns,1), [],[],'sp', 0.01, 'mar', 0.05);
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
    xlabel(h.ax(1), 'Age (months)');
    lbl = {...
            'w = von Bertalanffy body weight'
            'l = survivorship'
            'w*l'
            '$\sum_{0}^{a}{Z_a} - a\frac{BA}{B}$'};
    legendflex(h.ax(1), lbl, 'ref', h.ax(1), 'nrow', 1, 'anchor', {'n','s'}, ...
        'xscale', 0.5, 'interpreter', 'Latex', 'buffer', [0 5]);
end
    

%------------
% Main calcs
%------------

function [Out, D] = editstanzacalcs(a, k, bab, blead, qblead, z, da)
%
% a:        vector, age (in months) at the start of stanza
% k:        vonBertalanffy growth curve constant (annual)
% bab:      relative biomass accumulation rate (BA/B), (annual)
% blead:    biomass of leading stanza
% qblead:   consumption rate of leading stanza (/yr)
% da:       discretization interval (months)

% Setup of discretization.  Note that Ecopath uses 90% as the upper
% bound, but that seems to leave out a good amount, so I'm upping it to
% 99.99%.
    
amax = log(1 - 0.9999^(1/3))./-(k/12); 
amax90 = log(1 - 0.9^(1/3))./-(k/12); % Old one, used for plots
xa = 0:da:ceil(amax);

% Which stanza does each month-value fall into?

ia = sum(bsxfun(@gt, xa, a));
ia(ia == 0) = 1;

% BAB value, either scalar for original EwE formulation (one rate for
% whole group), or mutliple values for Rpath extension (one per age
% group)

if isscalar(bab)
    bab2 = bab./12; % to monthly
else
    bab2 = bab(ia)./12; % expand, and to monthly
end
    
% Expand Z values to month-ages

z = z(ia);

% Biomass and consumption of each stanza

zsum = cumsum(z./12 .* da);    % Integral of Z, 0 to a
l = exp(-zsum- xa'.*bab2);      % survivorship

num = l./sum(l);
w = (1 - exp(-k.*xa'./12)).^3; % weight    
wwa = w.^(2/3);

bsa = l.*w./(sum(l.*w));       % relative biomass at a given age
csa = (l.*wwa./(sum(l.*w)));   % relative consumption

alim = [a; Inf];
[bs, qs] = deal(zeros(size(a)));
for ia = 1:length(a)
    isin = xa >= alim(ia) & xa < alim(ia+1);
    bs(ia) = sum(l(isin).*w(isin))./sum(l.*w);

    qs(ia) = sum(l(isin).*wwa(isin))./sum(l.*wwa);

end

btot = blead./bs(end);
Out.b = btot .* bs;

% Note: I'm not getting the exact same values here as in EwE6. I
% think this might be related to the K calculation in the original
% code, but I haven't yet figured out exactly what they're doing
% here...
%
% (from Sub CalculateStanzaParameters, in cEcosimModel.vb)
%
%-----
% K = 0   'temporarily use k to sure the sum:
% For Age = first(BaseCB) To Second(BaseCB)
%     K = K + m_stanza.SplitNo(isp, Age) * m_stanza.WWa(isp, Age)
% Next
% 
% 
% If K > 0 Then K = cb(BaseCB) * Bio(BaseCB) / K 'THIS IS THE REAL CONSTANT k
%----
%
% I think this is just checking the the sum over what I call qs is 1,
% and normalizing if not.  So far my tests are always at 1 even without
% this check.  

qtot = blead.*qblead./qs(end);    
q = qtot .* qs;
Out.qb = q./Out.b;

% Calculate BA

Out.ba = Out.b.*bab;

% Values for plots

isless = xa < amax90;
D.x = xa(isless);
D.y = [w l w.*l zsum-xa'.*bab];
D.y = D.y(isless,:);





