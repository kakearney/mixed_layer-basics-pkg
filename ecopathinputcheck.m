function A = ecopathinputcheck(A, warnoff)
%ECOPATHINPUTCHECK Checks and fixes input in an Ewe input structure
% 
% B = ecopathinputcheck(A)
% B = ecopathinputcheck(A, warnoff)
%
% This function checks an Ewe input structure for proper dimensions and to
% verify certain values.  If incorrect values are found (for example,
% non-zero values for a primary producer's consumption/biomass ratio), the
% values are corrected.
%
% A note on multi-stanza groups: As noted in calcstanza.m, I do not always
% exactly reproduce the values of B and QB seen in EwE6.  For this reason,
% this function checks that multi-stanza group values are within 0.05% of
% the values calculated by my algorithm.  If B or QB data are missing for
% non-leading stanzas, that data is filled in.  If a stanza group has
% a biomass value outside the 0.05% error range, the biomass values for all
% stanzas of that particular group are recalculated and replaced; likewise
% for QB values.  Otherwise, they are left alone.
%
% Input variables:
%
%   A:          Ewe input structure
%
%   warnoff:    logical scalar.  If true, no warnings are issued when
%               corrections are made to the structure.  Default: false.
%
% Output variables:
%
%   B:          structure identical to A but with corrections made if
%               necessary.

% Copyright 2008-2015 Kelly Kearney

if nargin < 2
    warnoff = false;
end

if ~isscalar(A) || ~isstruct(A)
    error('Input structure must be scalar');
end

%----------------------------
% Check sizes
%----------------------------

group0 = {'ngroup', 'ngear', 'nlive'};
group1 = {'areafrac', 'b', 'pb', 'qb', 'ee', 'ge', 'gs', 'dtImp', 'bh', ...
          'pp', 'immig', 'emig', 'emigRate', 'ba', 'baRate'};

tf = isfield(A, [group0 group1]);
if ~all(tf)
    tmp = [group0 group1];
    str = sprintf('%s, ', tmp{~tf});
    error('Missing fields in input: %s', str(1:end-2));
end
      
if ~all(cellfun(@(x) isscalar(x) && isnumeric(x), {A.ngroup, A.ngear, A.nlive}))
    error('ngroup, nlive, and ngear must be numeric scalars');
end

group1size = cellfun(@(x) size(A.(x)), group1, 'uni', 0);

if ~isequal([A.ngroup 1], group1size{:})
    errlist = sprintf('  %s, %s, %s, %s, %s, %s, %s, %s,\n', group1{:});
    errlist(end-1:end) = [];
    error('The following fields must be ngroup x 1 vectors: \n%s', errlist);
end

if ~isequal([A.ngroup A.ngroup-A.nlive], size(A.df))
    error('df must be ngroup x ndet array');
end

if ~isequal([A.ngroup A.ngroup], size(A.dc))
    error('dc must be ngroup x ngroup array');
end

if ~isequal([A.ngroup A.ngear], size(A.landing), size(A.discard))
    error('landing and discard must be ngroup x ngear arrays');
end

if ~isequal([A.ngear A.ngroup-A.nlive], size(A.discardFate))
    error('discardFate must be ngear x ndet array');
end

%----------------------------
% Check content
%----------------------------

hasdtimp = (1:A.ngroup)' <= A.nlive & A.dtImp > 0;
if any(hasdtimp)
    if ~warnoff
        warnmessage('Non-zero values found for detritus import of live groups;\n values have been reset to 0', ...
                hasdtimp, 'dtImp', A);
    end
    A.dtImp(hasdtimp) = 0;
end

nandtimp = (1:A.ngroup)' > A.nlive & isnan(A.dtImp);
if any(nandtimp)
    if ~warnoff
        warnmessage('NaN found in detritus import, replacing with 0', ...
                nandtimp, 'dtImp', A);
    end
    A.dtImp(nandtimp) = 0;
end

epfields = {'b', 'pb', 'qb', 'ee', 'ge', 'bh'};
for ifield = 1:length(epfields)
    isneg = (A.(epfields{ifield}) < 0);
    if any(isneg)
        if ~warnoff
            msg = sprintf('Negative placeholders found in %s field, replacing with NaN', epfields{ifield});
            warnmessage(msg, ...
                isneg, epfields{ifield}, A);
        end
        A.(epfields{ifield})(isneg) = NaN;
    end
end
    
qbnotzero = A.pp >= 1 & A.qb ~= 0;
if any(qbnotzero)
    if ~warnoff
        warnmessage('Non-zero value found for a producer or detrital Q/B, replacing with zero', ...
                qbnotzero, 'qb', A);
    end
    A.qb(qbnotzero) = 0;
end

genotzero = A.pp >= 1 & A.ge ~= 0;
if any(genotzero)
    if ~warnoff
        warnmessage('Non-zero value found for a producer or detrital GE, replacing with zero', ...
                genotzero, 'ge', A);
    end
    A.ge(genotzero) = 0;
end

gsnotzero = A.pp >= 1 & A.gs ~= 0;
if any(gsnotzero)
    if ~warnoff
        warnmessage('Non-zero value found for a producer or detrital GS, replacing with zero', ...
                gsnotzero, 'gs', A);
    end
    A.gs(gsnotzero) = 0;
end

pbnotzero = A.pp == 2 & A.pb ~= 0;
if any(pbnotzero)
    if ~warnoff
        warnmessage('Non-zero value found for detrital P/B, replacing with zero', ...
                pbnotzero, 'pb', A);
    end
    A.pb(pbnotzero) = 0;
end

detbmissing = isnan(A.b) & A.pp == 2;
if any(detbmissing)
    if ~warnoff
        warnmessage('Detritus groups found with missing biomass; replacing with 0', ...
                detbmissing, 'b', A);
    end
    A.b(detbmissing) = 0;
end

% Check diet: should sum to 1 for non-producers.  Any missing goes to
% import.

% tol = 1e-6;
% dcsum = sum(A.dc,1);
% dctarget = max(1 - A.pp,0);
% 
% ppadjust = (1-dcsum')<tol & (A.pp < 1 & A.pp > 0);
% A.dc(:,ppadjust) = bsxfun(@times, A.dc(:,ppadjust), dctarget(ppadjust)');

% Added input field in later versions of code, derive from DC if not
% present

if ~isfield(A, 'import')
    A.import = 1 - sum(A.dc,1)' - min(A.pp,1);
end

%----------------------------
% Check Ecosim-related 
% content (optional fields)
%----------------------------

groupinfo = {'maxrelpb', 'maxrelfeed', 'feedadj', 'fracsens', ...
            'predeffect', 'densecatch', 'qbmaxqb0', 'switchpower'};

% Check group info sizes.  Group info fields are listed in GUI without
% detritus, so input coming from these tables may be too short. 

for ifd = 1:length(groupinfo)
    if isfield(A, groupinfo{ifd})
        if isequal(size(A.(groupinfo{ifd})), [A.ngroup-1 1])
            A.(groupinfo{ifd}) = [A.(groupinfo{ifd}); zeros(A.ngroup-A.nlive,1)];
        elseif ~isequal(size(A.(groupinfo{ifd})), [A.ngroup 1])
            error('%s field must be ngroup x 1 vector', groupinfo{ifd});
        end
    end
end
  
% Max rel. P/B should be 0 for non-producers
    
if isfield(A, 'maxrelpb')
    maxpbnotzero = A.pp ~= 1 & A.maxrelpb ~= 0;
    if any(maxpbnotzero)
        if ~warnoff
            warnmessage('Non-zero max rel P/B found for non-producers, replacing with 0', ...
                maxpbnotzero, 'maxrelpb', A);
        end
        A.maxrelpb(maxpbnotzero) = 0;
    end
end

% Other group info should be 0 for producers and detritus

for ifd = 2:length(groupinfo)
    if isfield(A, groupinfo{ifd})
        ginotzero = A.pp ~= 0 & A.(groupinfo{ifd}) ~= 0;
        if any(ginotzero)
            if ~warnoff
                msg = sprintf('Non-zero %s found for producer or detritus, replacing with 0', groupinfo{ifd});
                warnmessage(msg, ...
                    ginotzero, groupinfo{ifd}, A);
            end
            A.(groupinfo{ifd})(ginotzero) = 0;
        end
    end
end

% kv should be 0 where dc is 0

if isfield(A, 'kv')
    if ~isequal(size(A.kv), [A.ngroup A.ngroup])
        error('kv field mst be ngroup x ngroup array');
    end
    kvnotzero = A.dc == 0 & A.kv ~= 0;
    if any(kvnotzero(:))
        if ~warnoff
            warnmessage('Non-zero kv found for non-feeding links, replacing with 0', ...
                kvnotzero, 'kv', A);
        end
        A.kv(kvnotzero) = 0;
    end
end
    
%----------------------------
% Check multi-stanza values
%----------------------------

% Check for appropriate variables.

stanzavars = {'stanza', 'vbK', 'stanzadata'};
hasstanzas = all(isfield(A, stanzavars)) && any(A.stanza > 0);

% Check that multi-stanza group values are consistent with each other if
% filled in already.  Fill in if not.

if hasstanzas
    Tmp = calcstanza(A);
    
    % If non-leading stanza group data was missing, fill it in (all or
    % nothing... if one stanza-group is missing, all other stanzas of that
    % group need to be recalculated)
    
    bfill = isnan(A.b) & ~isnan(Tmp.b);
    bchange = ismember(A.stanza, unique(A.stanza(bfill)));
    A.b(bchange) = Tmp.b(bchange);
    qfill = isnan(A.qb) & ~isnan(Tmp.qb);
    qchange = ismember(A.stanza, unique(A.stanza(qfill)));
    A.qb(qchange) = Tmp.qb(qchange);
    
    % Check for any other changes.  Keep original data if it's within my
    % tolerance (meaning probably correct, just picking up the differences
    % between my implementation of the staza calculations vs EwE6's
    % implementation).  If it's further off, assume incorrect data, and
    % replace data for all stanzas of that group.
    
    changeflag = iscell(A.stanzadata.BABsplit) || ~any(isnan(A.stanzadata.BABsplit));
    
    berr = (Tmp.b - A.b)./A.b;
    qerr = (Tmp.qb - A.qb)./A.qb;
    baerr = (Tmp.ba - A.ba)./A.ba;
    
    tol = 0.005;
    bwrong = abs(berr) > tol;
    qwrong = abs(qerr) > tol;
    bawrong = abs(baerr) > tol;
    
    if changeflag
        if any([bwrong; qwrong; bawrong])

            bchange = ismember(A.stanza, unique(A.stanza(bwrong)));

            if ~warnoff
                msg = 'Multi-stanza data inconsistent; replacing all B, Q/B, and BA in stanza groups related to these B...';
                tmp = [A.name(bwrong) num2cell([A.b(bwrong) Tmp.b(bwrong)])]';
                str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                warning('%s:\n%s', msg, str);
            end
            A.b(bchange) = Tmp.b(bchange);

            qchange = ismember(A.stanza, unique(A.stanza(qwrong)));

            if ~warnoff
                msg = '...Q/B...';
                tmp = [A.name(qwrong) num2cell([A.qb(qwrong) Tmp.qb(qwrong)])]';
                str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                warning('%s:\n%s', msg, str);
            end
            A.qb(qchange) = Tmp.qb(qchange);
            
            bachange = ismember(A.stanza, unique(A.stanza(bawrong)));
            if ~warnoff
                msg = '... and BA';
                tmp = [A.name(bawrong) num2cell([A.ba(bawrong) Tmp.ba(bawrong)])]';
                str = sprintf('  %s: %.4g -> %.4g\n', tmp{:});
                warning('%s:\n%s', msg, str);
            end
            A.ba(bachange) = Tmp.ba(bachange);
            
        end   
    else
        if ~warnoff
            warning('Not checking stanza parameters');
        end
    end
end

% Warning messages

function warnmessage(msg, mask, fld, A)
tmp = [A.name(mask) num2cell(A.(fld)(mask))]';
str = sprintf('  %s (%.2f)\n', tmp{:});
warning('%s:\n%s', msg, str);

    
    
    

