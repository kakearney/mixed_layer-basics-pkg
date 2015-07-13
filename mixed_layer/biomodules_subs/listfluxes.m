function list = listfluxes(type, Idx, links)
%LISTFLUXES Return indices of fluxes in nemurokak or wce model
%
% list = listfluxes('nemurokak', Idx)
% list = listfluxes('wce', Idx, links)
%
% Input variables:
%
%   Idx:    1 x 1 structure with fields (ps, pl, zs, zl, etc) indicating
%           the index of the biological state variables corresponding to
%           each nemuro-derived variable (see wce and nemurokak setup)
%
%   links:  ng x ng array indicating type of predator/prey interaction
%           1 = zooplankton eat zooplankton
%           2 = nekton eat zooplankton
%           3 = nekton eat nekton
%           4 = zooplankton eat phytoplankton
%           Note: For mortality flux indices, this function assumes that
%           all living critters are involved in at least one predator-prey
%           interaction.  If not true, I'll need to update this.
%
% Output variables:
%
%   list:   n x 3 cell array, where column 1 indicates the type of flux,
%           column 2 source group index, and column 3 the sink group index

% Copyright 2011 Kelly Kearney

switch type
    
    case 'nemurokak'
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi
            Idx.fe      Idx.psfe
            Idx.fe      Idx.plfe];

        A.exx = [...
            Idx.ps      Idx.don
            Idx.pl      Idx.don
            Idx.plsi    Idx.sioh4];
%             Idx.mys     Idx.fe];

        A.res = [...
            Idx.ps      Idx.no3     
            Idx.ps      Idx.nh4     
            Idx.pl      Idx.no3     
            Idx.pl      Idx.nh4     
            Idx.plsi    Idx.sioh4];   
%             Idx.mys     Idx.fe];
        
        A.gra = [...
            Idx.ps      Idx.zs
            Idx.ps      Idx.zl
            Idx.pl      Idx.zl
            Idx.pl      Idx.zp
            Idx.zs      Idx.zl
            Idx.zs      Idx.zp
            Idx.zl      Idx.zp
            Idx.psfe    Idx.mys
            Idx.plfe    Idx.mys];
        
        A.pre = [];
        
        A.exc = [...
            Idx.zs      Idx.nh4
            Idx.zl      Idx.nh4
            Idx.zp      Idx.nh4
            Idx.psfe    Idx.mys
            Idx.plfe    Idx.mys];
%             Idx.mys     Idx.fe];
        
        A.ege = [...
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal
            Idx.psfe    Idx.pofe
            Idx.plfe    Idx.pofe];
%             Idx.mys     Idx.fe];
        
        A.mor = [...
            Idx.ps      Idx.pon
            Idx.pl      Idx.pon
            Idx.zs      Idx.pon
            Idx.zl      Idx.pon
            Idx.zp      Idx.pon
            Idx.plsi    Idx.opal
            Idx.psfe    Idx.pofe
            Idx.plfe    Idx.pofe];
%             Idx.mys     Idx.fe];
        
        A.dec = [...
            Idx.nh4     Idx.no3
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
            Idx.pon     Idx.don
            Idx.opal    Idx.sioh4
            Idx.pofe    Idx.fe
            Idx.fe      Idx.pofe];
        
    case 'wce'
        
        A.gpp = [...
            Idx.no3     Idx.ps
            Idx.nh4     Idx.ps
            Idx.no3     Idx.pl
            Idx.nh4     Idx.pl
            Idx.sioh4   Idx.plsi
            Idx.fe      Idx.psfe
            Idx.fe      Idx.plfe];
        
        A.npp = [...
            Idx.mys     Idx.ps
            Idx.mys     Idx.pl];

        A.exx = [...
            Idx.ps      Idx.don
            Idx.pl      Idx.don
            Idx.plsi    Idx.sioh4];
%             Idx.mys     Idx.fe];

        A.res = [...
            Idx.ps      Idx.no3     
            Idx.ps      Idx.nh4     
            Idx.pl      Idx.no3     
            Idx.pl      Idx.nh4     
            Idx.plsi    Idx.sioh4];   
%             Idx.mys     Idx.fe];
        
        [pry1, prd1] = find(links == 2 | links == 3); % nektonic
        [pry2, prd2] = find(links == 1 | links == 4); % planktonic
        
        A.gra = [...
            pry2        prd2];
        
        A.pre = [...
            pry1        prd1];
        
        
        unqprd = unique([prd1; prd2]);
        
        A.exc = [...
            unqprd      Idx.nh4.*ones(size(unqprd))
            Idx.psfe    Idx.mys
            Idx.plfe    Idx.mys];
%             Idx.mys     Idx.fe];
        
        A.ege = [...
            unqprd      Idx.pon.*ones(size(unqprd))
            Idx.plsi    Idx.opal
            Idx.psfe    Idx.pofe
            Idx.plfe    Idx.pofe];
%             Idx.mys     Idx.fe];
        
        alllive = unique([pry1; pry2; prd1; prd2]); % Assumes all live guys eat or are eaten
        
        A.mor = [...
            alllive     Idx.pon.*ones(size(alllive))
            Idx.plsi    Idx.opal
            Idx.psfe    Idx.pofe
            Idx.plfe    Idx.pofe];
%             Idx.mys     Idx.fe];
        
        A.dec = [...
            Idx.nh4     Idx.no3
            Idx.don     Idx.nh4
            Idx.pon     Idx.nh4
            Idx.pon     Idx.don
            Idx.opal    Idx.sioh4
            Idx.pofe    Idx.fe
            Idx.fe      Idx.pofe];
        
    otherwise
        error('Unrecognized option for biological module');
        
        
end
        
% Reformat into list

list = cell(0,3);

fld = fieldnames(A);
for ifld = 1:length(fld)
    nlnk = size(A.(fld{ifld}),1);
    
    list = [list; repmat(fld(ifld), nlnk, 1), num2cell(A.(fld{ifld}))];
end
            
        
