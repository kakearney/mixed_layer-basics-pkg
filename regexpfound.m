function found = regexpfound(str, expression);
%REGEXPFOUND Determine whether a regular expression is in a string
%
% found = regexpfound(str, expression);
% found = regexpfound(cells, expression);
%
% This function returns a logical array indicating whether the expression
% was found in the input string or cell array.  It acts similar to regexp
% but with output of a logical array rather than the typical cell array.
%
% Input variables:
%
%   str:        vector character array
%
%   cells:      cell array of strings
%
%   expression: regular expression to match
%
% Output variables:
%
%   found:      logical array same size as cells (1 x 1 if input was a
%               string), true if expression was found 1 or more times in
%               the string

% Copyright 2008 Kelly Kearney

if ischar(str)
    str = cellstr2(str);
end

found = cellfun(@(x) ~isempty(x), regexp(str, expression));