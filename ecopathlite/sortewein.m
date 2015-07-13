function [Ewein, isrt] = sortewein(Ewein)
%SORTEWEIN Sort groups in Ewe input structure by trophic level
%
% [New, isrt] = sortewein(Ewein)
%
% This function rearranges the groups in a Ewe input structure so that they
% are in descending order by trophic level (and also makes sure detritus
% groups are listed last).
%
% Input variables:
%
%   Ewein:  Ewe input structure
%
% Output variables
%
%   New:    Ewein input structure with groups sorted by trophic level
%
%   isrt:   index vector showing old positions of new groups, such that
%           old(isrt) = new

% Copyright 2009 Kelly Kearney


Ewein = ecopathinputcheck(Ewein, true);
if ~isfield(Ewein, 'name')
    nameflag = true;
    Ewein.name = cellstr(num2str((1:Ewein.ngroup)', 'Group %d'));
end
names = Ewein.name;

% First, check to make sure detritus groups come last.  If not, the trophic
% level calculations will be wrong (Note: I assume the detritus won't need
% to be reordered... is that always valid?)

[srt, isrtpp] = sort(Ewein.pp);
Ewein = resortewe(Ewein, isrtpp);

% Now sort by trophic level

tl = trophiclevel(Ewein.dc, Ewein.pp, Ewein.nlive, Ewein.ngroup);
[tlsrt, isrttl] = sort(tl, 'descend');

Ewein = resortewe(Ewein, isrttl);

% Match indices

[tf, isrt] = ismember(Ewein.name, names);

%----------------------
% Reorder Ewein data
%----------------------

function A = resortewe(A, isrt)

group1 = {'areafrac', 'b', 'pb', 'qb', 'ee', 'ge', 'gs', 'dtImp', 'bh', ...
          'pp', 'immig', 'emig', 'emigRate', 'ba', 'baRate', 'df', ...
          'landing', 'discard', 'maxrelpb', 'qbmaxqb0', 'name'};
group2 = {'dc', 'kv'};

for ig = 1:length(group1)
    if isfield(A, group1{ig})
        A.(group1{ig}) = A.(group1{ig})(isrt,:);
    end
end

for ig = 1:length(group2)
    if isfield(A, group2{ig})
        A.(group2{ig}) = A.(group2{ig})(isrt,isrt);
    end
end

