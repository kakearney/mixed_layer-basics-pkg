function tf = checkbalance(In)
%CHECKBALANCE Check whether Ecopath model is balanced
%
% tf = checkbalance(In)
%
% Output variables:
%
%   tf: logical scalar, true if balanced

In = ecopathinputcheck(In, true);
Ep = ecopathlite(In);

tf = all(Ep.ee <= 1 & Ep.ee >= 0);