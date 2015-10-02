function A = calcstanza(A)
% CALCSTANZA Calculate B and Q/B values for multi-stanza Ecopath groups
%
% A = calcstanza(In)
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
%   In: Ecopath input structure.  B and QB of non-leading multi-stanza 
%       groups can be unknown (i.e. NaN).  Stanza-related fields (stanza,
%       stanzadata.stanzaID, stanzadata.BABsplit, ageStart, and vbK) must
%       be present. 
%
% Output variables:
%
%   A:  Ecopath input structure, identical to In except that non-leading
%       multi-stanza group B and QB values have been recalculated/filled
%       in.

% Copyright 2014-2015 Kelly Kearney

ns = length(A.stanzadata.StanzaID);

for is = 1:ns
    
    idx = find(A.stanza == is);
    
    [a, isrt] = sort(A.ageStart(idx));
    idx = idx(isrt);
   
    % Parameters
    
    k = A.vbK(idx(1));               % Curvature parameter
    bab = A.stanzadata.BABsplit(is)./12; % biomass acc. per unit biomass, monthly
    
    % Setup of discretization.  Note that Ecopath uses 90% as the upper
    % bound, but that seems to leave out a good amount, so I'm upping it to
    % 99.99%.
    
    da = 1;
    amax = log(1 - 0.9999^(1/3))./-(k/12); 
    xa = 0:da:ceil(amax);
    
    % Biomass and consumption of each stanza
    
    ia = sum(bsxfun(@gt, xa, a));
    ia(ia == 0) = 1;
    z = A.pb(idx(ia));
   
    zsum = cumsum(z./12 .* da);    % Integral of Z, 0 to a
    l = exp(-zsum- xa'.*bab);      % survivorship
    
    num = l./sum(l);
    w = (1 - exp(-k.*xa'./12)).^3; % weight    
    wwa = w.^(2/3);
    
    bsa = l.*w./(sum(l.*w)); % relative biomass at a given age
    csa = (l.*wwa./(sum(l.*w))); % relative consumption
    
    alim = [a; Inf];
    [bs, qs] = deal(zeros(size(a)));
    for ia = 1:length(a)
        isin = xa >= alim(ia) & xa < alim(ia+1);
        bs(ia) = sum(l(isin).*w(isin))./sum(l.*w);
        
        qs(ia) = sum(l(isin).*wwa(isin))./sum(l.*wwa);
        
    end
    
    btot = A.b(idx(end))./bs(end);
    b = btot .* bs;
    A.b(idx) = b;
    
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

    qtot = A.b(idx(end)).*A.qb(idx(end))./qs(end);    
    q = qtot .* qs;
    A.qb(idx) = q./b;

    % Plot to check
    
    plotflag = false;
    if plotflag
        figure;

        [hl, ha] = plots(xa', [w l w.*l]);
        gridxy(a, [],'parent', ha(1));
        cmap = get(0, 'DefaultAxesColorOrder');
        col = num2cell(cmap(1:length(hl),:),2);
        set(hl, {'color'}, col);
        set(ha, {'ycolor'}, col);
        xlabel(ha(1), 'Age (months)');
        ylabel(ha(1), 'w = von Bertalanffy body weight')
        ylabel(ha(2), 'l = survivorship')
        ylabel(ha(3), 'w*l')
    end
    
end
