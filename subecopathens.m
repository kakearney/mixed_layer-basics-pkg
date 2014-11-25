function [NewIn, x] = subecopathens(Ewein, x, idx)
%SUBECOPATHENS Create Ecopath input structure from ensemble data
%
% NewIn = subecopathens(Ewein, x, idx)
%
% This function substitutes the results of createensemble back into an
% Ecopath input structure.  When substituting DC, it also renormalizes
% predator diet back to 1.  When substituting B and QB of multi-stanza
% groups, it only substitutes the values for the leading stanza, and
% calculates values for the remaining stanza groups based on this.
%
% Input variables:
%
%   Ewein:  original Ecopath input structure
%
%   x:      1 x 6 cell array, ensemble parameter data (see
%           createensemble.m) 
%
%   idx:    1 x 6 cell array, indices corresponding to parameter data (see
%           createensemble.m)
%
% Output variables:
%
%   NewIn:  nx x 1 Ecopath input structure, where nx corresponds to the
%           number of parameter sets in x (i.e. size(x{ii},1) ).
%
%   x:      same as input, except x{4} data (i.e. DC data) has been
%           normalized so predator diet sums to 1

% Copyright 2014 Kelly Kearney


nx = size(x{1},1);
ncol = cellfun(@length, idx);

hasstanza = isfield(Ewein, 'stanza') && any(Ewein.stanza > 0);

NewIn = repmat(Ewein, nx, 1);
for ii = 1:nx
    
    % Substitute values
    
    NewIn(ii).b( idx{1}) = x{1}(ii,:);
    NewIn(ii).pb(idx{2}) = x{2}(ii,:);
    NewIn(ii).qb(idx{3}) = x{3}(ii,:);
    NewIn(ii).dc(idx{4}) = x{4}(ii,:);
    NewIn(ii).ee(idx{5}) = x{5}(ii,:);
    NewIn(ii).ge(idx{6}) = x{6}(ii,:);
    
    NewIn(ii).bh = NewIn(ii).b ./ NewIn(ii).areafrac;
    
    % Adjust multi-stanza groups if necessary
    
    if hasstanza
        NewIn(ii) = calcstanza(NewIn(ii));
    end
    
    % Make sure DC sums to 1 for all predators
    
    NewIn(ii).dc = bsxfun(@rdivide, NewIn(ii).dc, sum(NewIn(ii).dc,1));
    NewIn(ii).dc(isnan(NewIn(ii).dc)) = 0;
    
    x{4}(ii,:) = NewIn(ii).dc(idx{4});
     
end


    