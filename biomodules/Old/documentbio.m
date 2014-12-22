function documentbio
%DOCUMENTBIO
%
% This function reformats the info in biomodules.txt into the help format
% used by all biological modules.  Easier than cutting, pasting, typing,
% and formatting every time.  Text saved to biomodulehelp.txt.

fid = fopen('biomodule.m');
bmtxt = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
bmtxt = bmtxt{1};
idx1 = find(strncmp('% Input variables for ''init'' mode', bmtxt, 33));
idx2 = find(strncmp('% Copyright', bmtxt, 11));
idx3 = find(regexpfound(bmtxt, '!!! LIST USER INPUT VARIABLES HERE !!!'));

maintxt1 = bmtxt(idx1:idx3-2);
maintxt2 = bmtxt(idx3+1:idx2);
indent = strfind(bmtxt{idx3}, '!');
indent = indent(1);

% Read biomodule descriptions

fid = fopen('biomodules.txt');
txt = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
txt = txt{1};

% Reformat

idx = find(strcmp(txt, 'Name:'));
nmod = length(idx);
name = txt(idx+1);
h1 = txt(idx+3);
dsc1 = txt(idx+5);

descrip = txt(idx+5);
for ii = 1:nmod
    x1 = idx(ii) + 5;
    if ii == nmod
        x2 = length(txt);
    else
        x2 = idx(ii+1) - 1;
    end
    dscrAndVar = txt(x1:x2);
    varidx = find(strncmp('Vars:', dscrAndVar, 5));
    descrip{ii} = dscrAndVar(1:varidx-1);
    vartxt{ii} = dscrAndVar(varidx+1:end);

end

tabstart = 1:4:75;

fulltxt = []; 
for im = 1:nmod
    
    h1txt = sprintf('%%%s %s', upper(name{im}), h1{im});
    descriptxt = wraptochar(descrip{im}, 73, 2);
    usetxt = createusetxt(name{im});
    
    isemp = cellfun(@isempty, vartxt{im});
    newvar = find(~isemp & [true; isemp(1:end-1)]);
    vars = vartxt{im}(newvar);
    nchar = max(cellfun(@length, vars));
    col1 = tabstart(find(tabstart > (nchar + indent), 1, 'first'));
    varlist = [];
    for iv = 1:length(vars)
        if iv ~= length(vars)
            vardef = vartxt{im}(newvar(iv)+1:newvar(iv+1)-1);
        else
            vardef = vartxt{im}(newvar(iv)+1:end-1);
        end
        vardef = wraptochar(vardef, 75 - col1 + 1, col1-1);
        texttoinsert = [vars{iv} ':'];
        
        vardef{1}(indent:(indent+length(vars{iv}))) = [vars{iv} ':'];
        
        varlist = [varlist; {'%'}; vardef];
       
    end
    
    modtxt = [...
        h1txt
        '%'
        usetxt
        '%'
        descriptxt
        '%'
        maintxt1
        varlist
        maintxt2];
    
    fulltxt = [fulltxt; {''}; {'********************'}; {''}; modtxt];

end

printtextarray(fulltxt, 'biomodulehelp.txt');
    
function str = createusetxt(name)
str = {...
'% [bioinit, ismixed, bottomval, Biovars, names] = ...'
sprintf('%%    %s(''init'', In, Grd);', name)
'% [newbio, diag] = ...'
sprintf('%%    %s(''sourcesink'', oldbio, meanqi, temp, z, dz, Biovars, t, dt);', name)
'% wsink = ...'
sprintf('%%    %s(''vertmove'', oldbio, meanqi, temp, z, dz, Biovars, t, dt);', name)};


function b = wraptochar(str, nchar, padleft)
if ischar(str)
    str = {str};
end
b = [];
for il = 1:length(str)
    a = str{il};
    while length(a) > 0
        sp = [regexp(a, '\s') length(a)+1];
        idx = sp(find(sp < nchar, 1, 'last'));
        b = [b; {a(1:idx-1)}];
        a = strtrim(a(idx:end));
    end
end
for il = 1:length(b)
    b{il} = ['%' repmat(' ', 1, padleft-1) b{il}];
end
    
