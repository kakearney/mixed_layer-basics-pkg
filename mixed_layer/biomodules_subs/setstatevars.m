function [names, nbsv, nemidx, A] = setstatevars(BioIn)
%SETSTATEVARS State variable setup for mixed_layer wce/nemurokak modules

%-------------------------
% Set state variable names
%-------------------------

nemnames = {...
    'PS',       'Small Phytoplankton',          'molN/m^3'
    'PL',       'Large Phytoplankton',          'molN/m^3'
    'ZS',       'Small Zooplankton',            'molN/m^3'
    'ZL',       'Large Zooplankton',            'molN/m^3'
    'ZP',       'Pradatory Zooplankton',        'molN/m^3'
    'NO3',      'Nitrate',                      'molN/m^3'
    'NH4',      'Ammmonium',                    'molN/m^3'
    'PON',      'Particulate Organic Nitrogen', 'molN/m^3'
    'DON',      'dissolved Organic Nitrogen',   'molN/m^3'
    'SiOH4',    'Silicate',                     'molSi/m^3'
    'Opal',     'Particulate Opal',             'molSi/m^3'
    'Fe'        'Dissolved Iron'                'molFe/m^3'
    'PLsi'      'Large Phytoplankton Silica',   'molSi/m^3'
    'PSfe'      'Small Phytoplankton Iron'      'molFe/m^3'
    'PLfe'      'Large Phytoplankton Iron'      'molFe/m^3'
    'POFe'      'Particulate iron'              'molFe/m^3'};

if BioIn.isnem
    names = nemnames;
    nbsv = size(names,1);
else
    % Figure out which of the Ecopath model groups correspond to these
    % nemuro-derived variables.

    [tf, loc] = ismember(lower(nemnames(:,1)), BioIn.types); 

    nbsv = BioIn.EM.ngroup + sum(~tf);

    names = cell(nbsv,3);
    
    names(loc(tf),:) = nemnames(tf,:);
    names(loc(tf),2) = BioIn.EM.name(loc(tf));
    names(BioIn.EM.ngroup+1:end,:) = nemnames(~tf,:);

    % Name the remaining nekton and zooplankton groups, using N# and Z# for
    % short names and Ecopath names for long names

    isemp = cellfun('isempty', names(:,2));
    names(isemp,2) = BioIn.EM.name(isemp);

    [names{isemp,3}] = deal('mol N m^-3');
    isn = strcmp(BioIn.types, 'n');
    isz = strcmp(BioIn.types, 'z');

    names(isn,1) = cellstr(num2str((1:sum(isn))', 'N%02d'));
    names(isz,1) = cellstr(num2str((1:sum(isz))', 'Z%02d'));
end

%-------------------------
% Identifiers
%-------------------------

A.idx.ps    = find(strcmp(names(:,1), 'PS'));
A.idx.pl    = find(strcmp(names(:,1), 'PL'));
A.idx.no3   = find(strcmp(names(:,1), 'NO3'));
A.idx.nh4   = find(strcmp(names(:,1), 'NH4'));
A.idx.sioh4 = find(strcmp(names(:,1), 'SiOH4'));
A.idx.don   = find(strcmp(names(:,1), 'DON'));
A.idx.pon   = find(strcmp(names(:,1), 'PON'));
A.idx.opal  = find(strcmp(names(:,1), 'Opal'));
A.idx.fe    = find(strcmp(names(:,1), 'Fe'));
A.idx.plsi  = find(strcmp(names(:,1), 'PLsi'));
A.idx.psfe  = find(strcmp(names(:,1), 'PSfe'));
A.idx.plfe  = find(strcmp(names(:,1), 'PLfe'));
A.idx.zs    = find(strcmp(names(:,1), 'ZS'));
A.idx.zl    = find(strcmp(names(:,1), 'ZL'));
A.idx.zp    = find(strcmp(names(:,1), 'ZP'));
A.idx.pofe  = find(strcmp(names(:,1), 'POFe'));

%-------------------------
% Links (wce only)
%-------------------------

if ~BioIn.isnem
    
    % Marker for nemuro-derived variables, zoo, and nekton

    isnem = ismember(names(:,1), nemnames(:,1));
    [blah, nemidx] = ismember(nemnames(:,1), names(:,1));

    A.isextrazoo = regexpfound(names(:,1), 'Z[0-9][0-9]');
    A.isnek = regexpfound(names(:,1), 'N[0-9][0-9]');
    A.iszoo = ismember(names(:,1), {'ZS','ZL','ZP'}) | A.isextrazoo;
    A.isphy = ismember(names(:,1), {'PS','PL'});
    
    % Links
    
    A.links = zeros(nbsv);
    A.links(bsxfun(@and, A.iszoo, A.iszoo')) = 1; % Z-eat-Z
    A.links(bsxfun(@and, A.iszoo, A.isnek')) = 2; % N-eat-Z
    A.links(bsxfun(@and, A.isnek, A.isnek')) = 3; % N-eat-N
    A.links(bsxfun(@and, A.isphy, A.iszoo')) = 4; % Z-eat-P
    [r,c] = find(table2array(BioIn.EM.dc) == 0);
    idx = sub2ind(size(A.links), r, c);
    A.links(idx) = 0;
    
else
    nemidx = 1:nbsv;
end

%-------------------------
% Diapause-related
%-------------------------

A.diapause = BioIn.diapause;
if A.diapause
    names = [
        names
        {...
        'ZL1',       'Large Zooplankton (no diapause)', 'mol N m^-3'
        'ZL2',       'Large Zooplankton (diapause)',    'mol N m^-3'}];
    A.idx.zl1 = nbsv+1;
    A.idx.zl2 = nbsv+2;
    
    nbsv = nbsv + 2;
    
    if ~BioIn.isnem
        linktmp = zeros(nbsv);
        linktmp(1:nbsv-2,1:nbsv-2) = A.links;
        A.links = linktmp;

        A.isnek = [A.isnek; false; false];
    end

end

%-------------------------
% Non-state-variable 
% identifiers
%-------------------------

A.idx.mys   = nbsv + 1; % "Mystery" box, for source/sink terms coming from or going to nowhere
A.idx.fish  = nbsv + 2; % Fisheries 




    
