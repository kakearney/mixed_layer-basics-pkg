function displayecopath(varargin)
%DISPLAYECOPATH Prints ecopath results to screen
%
% displayecopath(In)
% displayecopath(In, Ep)
%
% This function displays the results of an ecopathlite simulation in the
% command window, using a table structure similar to the way Ecopath
% results are presented in the original EwE "Basic Estimates" gui.
%
% Note: I used to also try to replicate other panels (Key Indices,
% Mortalities, etc.) but I never really used those panels beyond the
% initial development phase, and the display relied on an outdated
% print-to-screen function, so I eliminated that behavior.   
%
% Input variables:
%
%   In:     Ewe input structure
% 
%   Ep:     ecopathlite results structure (see ecopathlite.m for details).
%           If not included, values set to NaN in the input structure will
%           be displayed as empty.

% Copyright 2012 Kelly Kearney

%---------------------------
% Parse input
%---------------------------

if nargin == 1
    epflag = false;
    In = varargin{1};
else
    epflag = true;
    In = varargin{1};
    Ep = varargin{2};
end


% 
% 
% 
% types = {...
%     'basic'     'Basic Estimates'
%     'key'       'Key Indices'
%     'mort'      'Mortalities'
%     'predmort'  'Predation Mortality'
%     'cons'      'Consumption'
%     'resp'      'Respiration'
%     'search'    'Search rates (Lotka-Volterra)'
%     'in'        'Input'};
% 
% needin = 8;
% needout = 1:7;
% 
% infield  = {'areafrac', 'b', 'pb', 'qb', 'ee', 'ge', 'ba', 'gs', 'dtImp'};
% outfield = {'ngroup', 'trophic', 'areafrac', 'bh', 'b', 'pb', 'qb', ...
%             'ee', 'ge', 'ba', 'baRate', 'migration', 'flowtodet', ...
%             'fishMortRate', 'predMortRate', 'migrationRate', ...
%             'otherMortRate', 'predMort', 'q0', 'q0Sum', 'respiration', ...
%             'searchRate', 'detexport'};
% 
% isin  = cellfun(@(x) isstruct(x) && all(isfield(x, infield)),  varargin);
% isout = cellfun(@(x) isstruct(x) && all(isfield(x, outfield)), varargin);
% istype = cellfun(@(x) ischar(x) && ismember(x, types(:,1)), varargin);
% 
% isdisp = ismember(types(:,1), varargin(istype));
% 
% if ~any(isdisp)
%     if any(isin)
%         isdisp(needin) = true;
%     end
%     if any(isout)
%         isdisp(needout) = true;
%     end
% end
% 
% if any(isdisp(needin))
%     if ~any(isin)
%         error('Must provide an Ewe input structure to display input');
%     end
%     A = varargin{isin};
% end
% 
% if any(isdisp(needout))
%     if ~any(isout)
%         error('Must provide an Ecopath output structure to display output');
%     end
%     S = varargin{isout};
% end
%     
% if exist('A', 'var') && isfield(A, 'name')
%     names = [A.name];
% else
%     if exist('S', 'var')
%         ng = S.ngroup;
%     else
%         ng = length(A.b);
%     end
%     names = [cellstr(num2str((1:ng)'))];
% end
% 
% %---------------------------
% % Collect data to display
% %---------------------------
% 
% % A display similar to 
% 
% if isdisp(1)
%     data{1} = [...
%     {'Trophic', 'Area', 'BH', 'B', 'P/B', 'Q/B', 'EE', 'GE'}
%     {'-','-', 'M/A', 'M/A', '/T', '/T', '-', '-'}
%     num2cell([S.trophic S.areafrac S.bh S.b S.pb S.qb S.ee S.ge])];
% end
% 
% if isdisp(2)
%     data{2} = [...
%         {'BA', 'BA rate', 'Migration', 'Det flow'}
%         {'M/A/T', '/T', 'M/A/T', 'M/A/T'}
%         num2cell([S.ba S.baRate S.migration S.flowtodet(1:S.ngroup)])
%         ];
% end
% 
% if isdisp(3)
%     data{3} = [...
%         {'Fishing', 'Predation', 'Migration', 'Other'}
%         num2cell([sum(S.fishMortRate,2) S.predMortRate S.migrationRate S.otherMortRate])
%         ];
% end
% 
% if isdisp(4)
%     data{4} = num2cell(S.predMort);
% end
% 
% if isdisp(5)
%     data{5} = num2cell([S.q0; sum(S.q0,1)]);
% end
% 
% if isdisp(6)
%     data{6} = num2cell(S.respiration);
% end
% 
% if isdisp(7)
%     data{7} = num2cell(S.searchRate);
% end
% 
% if isdisp(8)
%     data{8} = [...
%     {'Area', 'BH', 'B', 'P/B', 'Q/B', 'EE', 'GE'}
%     num2cell([A.areafrac A.bh A.b A.pb A.qb A.ee A.ge])];
% end
% 
% %---------------------------
% % Print to screen
% %---------------------------
% 
% for itype = find(isdisp')
%     
%     disp(' ');
%     disp(types{itype,2});
%     disp(' ');
%     
%     if (size(data{itype},1) - length(names)) == 0
%         disp(cprintf([names data{itype}], '-n', '%.4g', '-d', '|', '-N', ''));
%     else
%         namepad = cellstr(repmat(' ', size(data{itype},1) - length(names),1));
%         disp(cprintf([[namepad; names] data{itype}], '-n', '%.4g', '-d', '|', '-N', ''));
%     end
% end

% NEW

% Display basic input/estimates (TL HA BH B PB QB EE GE GS DI)

blank = nan(In.ngroup,1);

basicin = [blank ...
           In.areafrac ...
           blank ...
           In.b ...
           In.pb ...
           In.qb ...
           In.ee ...
           In.ge ...
           In.gs ...
           In.dtImp];
       
if epflag
    basicout = [Ep.trophic ...
               Ep.areafrac ...
               Ep.bh ...
               Ep.b ...
               Ep.pb ...
               Ep.qb ...
               Ep.ee ...
               Ep.ge ...
               blank ...
               blank];
else
    basicout = nan(In.ngroup, 10);
end

isfilled = isnan(basicin) & ~isnan(basicout);
isoutofrange = false(size(isfilled));
isoutofrange(:,7) = basicout(:,7) < 0 | basicout(:,7) > 1;
       
isinonly = isnan(basicout) & ~isnan(basicin);

dispval = basicin;
dispval(isfilled) = basicout(isfilled);

[nrow, ncol] = size(dispval);
    
gwidth = max(cellfun(@length, In.name));

dispstr = arrayfun(@(x) num2str(x,'%g'), dispval, 'uni', 0);
for ic = 1:ncol
    wdth(ic) = max(cellfun(@length, dispstr(:,ic)));
end
[dispstr{strcmp(dispstr, 'NaN')}] = deal(' ');

fmt = sprintf('%%-%ds ', 3, gwidth, wdth);
cprintf('-text', fmt, '   ', 'Name', 'TL', 'HA', 'BH', 'B', 'PB', 'QB', 'EE', 'GE', 'GS', 'DI');
fprintf('\n');

for ir = 1:nrow
    fmt = sprintf('%%2d: %%-%ds', gwidth+1);
    fprintf(fmt, ir, In.name{ir});
    for ic = 1:ncol
        fmt = sprintf('%%%ds', wdth(ic)+1);
        if isoutofrange(ir,ic)
            cprintf('red', fmt, dispstr{ir,ic});
        elseif isfilled(ir,ic)
            cprintf('blue', fmt, dispstr{ir,ic});
        else
            cprintf('text',fmt, dispstr{ir,ic});
        end
    end
    fprintf('\n');
end


           
           


