function [Out, D] = editstanzacalcs(a, k, bab, blead, qblead, z, da)
%EDITSTANZACALCS Replicate multi-stanza calculations from Ecopath
%
% [Out, D] = editstanzacalcs(a, k, bab, blead, qblead, z, da)
%
% This function replicates the EwE "Edit multi-stanza" calculations for a
% single species group.  It calculates B and QB for all mutli-stanza
% subgroups based on the B and QB of the leading stanza group.
%
% Input variables:
%
%   a:      n x 1 vector, where n is the number of subgroups (i.e. age
%           classes) in the multi-stanza group, holding age (in months) at
%           the start of stanza sub-group.  Values must be in order from
%           youngest to oldest stanza.
%
%   k:      scalar, vonBertalanffy growth curve constant (annual)
%
%   bab:    relative biomass accumulation rate (BA/B), (annual).  Can be
%           either a scalar, applying to all subgroups, or a n x 1 vector
%
%   blead:  scalar, biomass of leading stanza
%
%   qblead: scalar, consumption rate of leading stanza (/yr)
%
%   z:      n x 1 vector, mortality rate (/yr) of each subgroup, assumed to
%           be equal to the production rate (P/B) of the subgroup  
%
%   da:     discretization interval (months)
%
% Output variables
%
%   Out:    structure with the following fields:
%
%           b:  n x 1, biomass values for each subgroup
%
%           qb: n x 1, consumption rate (/yr) for each subgroup
%
%           ba: n x 1, biomass accumulation rate (/yr) for each subgroup
%
%   D:      structure holding plotting variables:
%
%           x:  nx x 1 vector, age, in months
%
%           y:  nx x 4 array.  Columns correspond to:
%               1:  w, individual weight at age, based on vonBertalanffy
%                   growth curve assuming weight propertional to length
%                   cubed.  
%               2:  l, survivorship at age, proportional to number at age
%               3:  w*l, population biomass at age
%               4:  exponent used in survivorship calculation, equal to
%                   sum_{0}^{a}{Z_a} - (a * (BA/B)_a), where a is age in
%                   months

% Copyright 2015 Kelly Kearney

% Setup of discretization.  Note that Ecopath uses 90% as the upper
% bound, and then has an accumulator function later on to deal with the 
% extra tail end of the biomass distribution.  I'm instead just using
% 99.99% as the upper limit, which gets me very close to the same numbers
% without adding computation (I don't loop over months like EwE6 does, so
% it's just a matter of having slightly larger matrices).
    
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
    bab = bab(:); % Make sure column vector
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
D.y = [w l w.*l zsum-xa'.*bab2];
D.y = D.y(isless,:);
