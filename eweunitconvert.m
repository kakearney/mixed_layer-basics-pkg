function varargout = eweunitconvert(oldunit, newunit, varargin)
%EWEUNITCONVERT Convert Ewe input structure between units
%
% [new1, new2, ...] = eweunitconvert(oldunit, newunit, old1, old2, ..., ...
%                                    param, val, ...)
%
% Input variables:
%
%   oldunit:    string specifying the units used in the old variable, in
%               the format M/A/T
%               M (mass):           'tons wet weight', 't ww', 'mmol N', 
%                                   'mmol n', 'mol N', 't C', 'g C'
%               A (area/volume):    'km^2', 'm^2', 'km^3', 'm^3'
%               T (time):           'year', 'yr', 'sec', 's'
%
%   newunit:    string specifying the units for the new variables, using
%               the same format as oldunit
%
%   old#:       variables to be converted.  These can be either Ewe input
%               structures, NemEp structures (see nemuroecopathdata) or
%               numeric arrays.  All numeric arrays are assumed to start
%               with the same unit type, specified by the 'varunit'
%               parameter
%
% Optional input variables (passed as parameter/value pairs), defaults in
% parentheses
%
%   wwCfrac:    ratio of carbon over wet weight (0.05)
%
%   c2n:        nitrogen/carbon ratio (16/106)
%
%   cmw:        molecular weight of carbon (12.0107)
%
%   depth:      depth of water column used for any area to volume
%               conversions, in km (NaN)
%
%   varunit:    type of unit used by input variables.  Can be either
%               'M/A/T' (flux rate), 'M/A' (concentration), or '1/T' (rate)
%               ('M/A/T')
%
% Output variables:
%
%   new#:       converted versions of all old# input variables

% Copyright 2008 Kelly Kearney

%------------------------------
% Parse input
%------------------------------

% Separate variables to be converted and parameter-value pairs

isparam = cellfun(@ischar, varargin);
if ~any(isparam)
    vars = varargin;
    pv = cell(0);
else
    pvstart = find(isparam, 1);
    vars = varargin(1:pvstart-1);
    pv = varargin(pvstart:end);
end

% Set optional parameters

Opt.wwCfrac = .05;                                      % C/wet weight for phytoplankton
Opt.c2n = 16/106; %17/133;                              % mol N/mol C
Opt.cmw = 12.0107;                                      % molecular weight of C
Opt.depth = NaN;                                        % depth for area to volume conversions, km
Opt.varunit = 'M/A/T';                                  % type of unit for non-input-output structures

Opt = parse_pv_pairs(Opt, pv);

% Parse unit strings

[m1,a1,t1] = strread(oldunit, '%s%s%s', 'delimiter', '/');
[m2,a2,t2] = strread(newunit, '%s%s%s', 'delimiter', '/');

%------------------------------
% Calculate conversion factors
%------------------------------

% Table of ratios

mass = {...
    'tons wet weight'   1
    't ww'              1
    't C'               Opt.wwCfrac
    'g C'               Opt.wwCfrac * 1e6
    'mmol N'            (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e9
    'mmol n'            (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e9
    'mol N'             (Opt.wwCfrac * Opt.c2n /Opt.cmw) * 1e6
    };

area = {...
    'km^2'              1
    'm^2'               1e6
    'km^3'              Opt.depth
    'm^3'               Opt.depth * 1e9
    };
    
time = {...
    'year'              1
    'yr'                1
    'day'               365
    'd'                 365
    'sec'               86400 * 365
    's'                 86400 * 365
    };

% Check input unit against list 

if ~all(ismember([m1 m2], mass(:,1)))
    error('Unrecognized mass unit');
end

if ~all(ismember([a1 a2], area(:,1)))
    error('Unrecognized mass unit');
end

if ~all(ismember([t1 t2], time(:,1)))
    error('Unrecognized mass unit');
end

% Old to new conversion factors

midx(1) = find(strcmp(m1, mass(:)));
midx(2) = find(strcmp(m2, mass(:)));

m2perm1 = mass{midx(2),2}/mass{midx(1),2};

aidx(1) = find(strcmp(a1, area(:)));
aidx(2) = find(strcmp(a2, area(:)));

a2pera1 = area{aidx(2),2}/area{aidx(1),2};

tidx(1) = find(strcmp(t1, time(:)));
tidx(2) = find(strcmp(t2, time(:)));

t2pert1 = time{tidx(2),2}/time{tidx(1),2};
 
% Field conversion factors

fac = {...
    '1/T'   1/t2pert1
    'M/A'   m2perm1/a2pera1
    'M/A/T' m2perm1/a2pera1/t2pert1
    };

%------------------------------
% Convert variables
%------------------------------

% Fields and units for Ecopath input structure

info = ecopathinputinfo;

% Fields and units for NEMURO-Ecopath structure

info2 = {...
    'b' 	'M/A'
    'qb'    '1/T'
    'q0'    'M/A/T'
    'q0sum' 'M/A/T'
    'pb'    '1/T'
    'ge'    'no unit'
    'gs'    'no unit'
    'dc'    'no unit'
    'm0'    'M/A/T'
    'm2'    'M/A/T'
    't'     'no unit'};
    

newvars = cell(size(vars));
for iv = 1:length(vars)
    
    isewein = false;
    isnem = false;
    if isstruct(vars{iv})
        flds = fieldnames(vars{iv});
        [tf, loc] = ismember(flds, info(:,1));
        
        if all(tf)
            isewein = true;
        else
            [tf, loc] = ismember(flds, info2(:,1));
            if all(tf)
                isnem = true;
            else
                error('Structure not recognized as either Ewein or NemEp');
            end
        end
    end
            
    
    if isewein % isstruct(vars{iv}) % Ewe input structure

%         flds = fieldnames(vars{iv});
%         
%         [tf, loc] = ismember(flds, info(:,1));
%         
%         if ~all(tf)
%             error('Input structure is not recognized as a EwE input structure');
%         end

        fldunit = info(loc, 3);

        New = vars{iv};
        for ifac = 1:length(fac)
            for ifld = 1:length(flds)
                if strcmp(fldunit{ifld}, fac{ifac,1})
                    New.(flds{ifld}) = vars{iv}.(flds{ifld}) .* fac{ifac,2};
                end
            end
        end
        
        newvars{iv} = New;
        
    elseif isnem
        
        fldunit = info2(loc,2);
        [blah, facidx] = ismember(fldunit, fac(:,1));
        
        New = vars{iv};
        for ifld = 1:length(flds)
            if facidx(ifld) ~= 0
                for ii = 1:numel(New.(flds{ifld}))
                    New.(flds{ifld}){ii} = vars{iv}.(flds{ifld}){ii} .* fac{facidx(ifld),2};
                end
            end
        end
        
        newvars{iv} = New;
        
    else % Other variables
        
        isfac = strcmp(Opt.varunit, fac(:,1));
        if ~any(isfac)
            error('Unrecognized variable unit format');
        end
        varfac = fac{isfac,2};
        newvars{iv} = vars{iv} .* varfac;
        
    end
end

varargout = newvars;
