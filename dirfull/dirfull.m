function varargout = dirfull(dirname)
%DIRFULL Same as dir but returns full filename with path
%
% dirfull(dirname)
% A = dirfull(dirname)
%
% Input variables:
%
%   dirname:    name of folder to list contents of
%
% Output variables:
%
%   A:          1 x nfile structure
%
%               name:       full path name of files
%
%               date:       modification date
%
%               bytes:      number of bytes allocated to the file
%
%               isdir:      1 if is a directory, 0 otherwise
%   
%               datenum:    modification date as a serial date number             

% Copyright 2013 Kelly Kearney

if nargout == 0
    dir(dirname);
else
    A = dir(dirname);

    dirpath = fileparts(dirname);

    for ia = 1:length(A)
        A(ia).name = fullfile(dirpath, A(ia).name);
    end

    varargout{1} = A;
end
