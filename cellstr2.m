function c = cellstr2(s)
%CELLSTR2 Create cell array of strings from character array.
%   C = CELLSTR2(S) places each row of the character array S into 
%   separate cells of C.
%
%   Use CHAR to convert back.
%
%   Another way to create a cell array of strings is by using the curly
%   braces: 
%      C = {'hello' 'yes' 'no' 'goodbye'};
%
%   This is a modified version of the cellstr function that does not
%   deblank the strings.  -Kelly Kearney
%
%   See also STRINGS, CHAR, ISCELLSTR.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.16.4.8 $  $Date: 2006/06/20 20:12:41 $
%==============================================================================

if ischar(s)
	if isempty(s)
	    c = {''};
	else
        if ndims(s)~=2
			error('MATLAB:cellstr:InputShape','S must be 2-D.')
        end
		[rows,cols]=size(s);%#ok ignore rows
		c = cell(rows,1);	
		for i=1:rows
            c{i} = s(i,:);
		end
	end
elseif iscellstr(s)
    c = s; 
else
	error('MATLAB:cellstr:InputClass','Input must be a string.')
end
%==============================================================================