function vars = nemvarnames(longstr)
%NEMVARNAMES Returns cell array holding NEMURO variable names
%
% vars = nemvarnames
% vars = nemvarnames('long')
% vars = nemvarnames('iron')
%
% Yes, I am that lazy
%
% Input variables:
%
%   'long': Return long names (e.g. 'Small Phytoplankton' as opposed to
%           'PS') 
%
%   'iron': Include the iron-related that I added to the original 11
%           variables.
%
% Output variables:
%
%   vars:   {'PS','PL','ZS','ZL','ZP','NO3','NH4','PON','DON','SiOH4',
%           'Opal'}, or a variation on that theme.
            

if nargin == 0
    longstr = 'original';
end

    
switch longstr

    case 'original'
        vars = {'PS','PL','ZS','ZL','ZP','NO3','NH4','PON','DON',...
                'SiOH4','Opal'};

    case 'long'
        vars = {'Small Phytoplankton', ...
                'Large Phytoplankton', ...
                'Small Zooplankton', ...
                'Large Zooplankton', ...
                'Predatory Zooplankton', ...
                'Nitrate', ...
                'Ammmonium', ...
                'Particulate Organic Nitrogen', ...
                'Dissolved Organic Nitrogen',  ...
                'Silicate', ...
                'Particulate Opal'};
    case 'iron'

        vars = {'PS','PL','ZS','ZL','ZP','NO3','NH4','PON','DON', ...
                'SiOH4','Opal', 'Fe', 'PLsi', 'PSfe', 'PLfe','POFe'};
            
end


