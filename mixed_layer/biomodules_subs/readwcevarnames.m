function [shortname, longname] = readwcevarnames(file)
%READWCEVARNAMES Read the names of biological state variables from file
%
% [shortname, longname] = readwcevarnames(file)
%
% Input variables:
%
%   file:   name of simulation output file
%
% Output vairables:
%
%   shortname:  short name of variables (used as field names when data is
%               read)
%
%   longname:   long name of variables (more descriptive)

% Copyright 2009 Kelly Kearney

[blah, blah, ext] = fileparts(file);
if isempty(ext)
    file = [file '.nc'];
end

Info = ncinfo(file);
varnames = {Info.Variables.Name};


% Old files

ng = sum(regexpfound(varnames, 'FG\d\d'));

if ng >= 1

    [shortname, longname] = deal(cell(ng,1));
    for ig = 1:ng
        shortname{ig} = sprintf('FG%02d', ig);
        longname{ig} = ncreadatt(file, shortname{ig}, 'long_name');
    end

else

% New files

    isvar = ismember(varnames, nemvarnames) | ...
            regexpfound(varnames, '^N\d\d$') | ...
            regexpfound(varnames, '^Z\d\d$') | ...
            ismember(varnames, {'Fe', 'PLsi', 'PSfe', 'PLfe', 'POFe'});
        
    shortname = varnames(isvar);
    longname = cell(size(shortname));
    for ig = 1:length(shortname)
        longname{ig} = ncreadatt(file, shortname{ig}, 'long_name');
    end
    
end

% close(nc);