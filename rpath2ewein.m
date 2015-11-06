function [A,B] = rpath2ewein(basename, varargin)
%RPATH2EWEIN Reads in Ecopath data from Rpath-formatted .csv files
%
% [A, B] = rpath2ewein(basename, p1, v1, ...)
%
% Input variables:
%
%   basename:   Base file name for the 4 input files
%
% Optional input variables (passed as paramete/value pairs)
%
%   basestr:    string appended to basename for basic data, not including
%               .csv extension ['_base']
%
%   dietstr:    string appended to basename for diet data, not including
%               .csv extension ['_diet']
%
%   juvsstr:    string appended to basename for stanza data, not including
%               .csv extension ['_juvs']
%
%   pedstr:     string appended to basename for pedigree data, not
%               including .csv extension ['_ped'] 
%
% Output variables:
%
%   A:          ecopath input structure (see ecopathlite.m)
%
%   B:          structure with full input tables from files

% Copyright 2015 Kelly Kearney

Opt.basestr = '_base';
Opt.dietstr = '_diet';
Opt.juvsstr = '_juvs';
Opt.pedstr = '_ped';

Opt = parsepv(Opt, varargin);

Base = readtable([basename Opt.basestr '.csv']);
Diet = readtable([basename Opt.dietstr '.csv'], 'ReadRowNames', true);
Juvs = readtable([basename Opt.juvsstr '.csv']);
Ped  = readtable([basename Opt.pedstr '.csv']);

B.Base = Base;
B.Diet = Diet;
B.Juvs = Juvs;
B.Ped = Ped;

bcol = parsecolname(Base);
dcol = parsecolname(Diet);
jcol = parsecolname(Juvs);
pcol = parsecolname(Ped);

isgroup = Base.Type <= 2;
isgear  = Base.Type == 3;
isdet   = Base.Type == 2;

% Basic

A.ngroup   = sum(isgroup);
A.nlive    = A.ngroup - sum(isdet);
A.ngear    = sum(isgear);
A.areafrac = ones(A.ngroup,1);
A.name     = strtrim(Base.Group(isgroup));
A.b        = Base.Biomass(isgroup);
A.pb       = Base.PB(isgroup);
A.qb       = Base.QB(isgroup);
A.ee       = Base.EE(isgroup);
A.ge       = Base.ProdCons(isgroup);
A.gs       = Base.Unassim(isgroup);
A.dtImp    = Base.DetInput(isgroup);
A.bh       = A.b;
A.pp       = Base.Type(isgroup);
A.ba       = Base.BioAcc(isgroup);
A.baRate   = zeros(A.ngroup,1);

A.immig    = zeros(A.ngroup,1);
A.emig     = zeros(A.ngroup,1);
A.emigRate = zeros(A.ngroup,1);

% Diet

A.dc = zeros(A.ngroup);

[tfr, locr] = ismember(Diet.Properties.RowNames, A.name);
[tfc, locc] = ismember(dcol, A.name);

A.dc(locr(tfr), locc(tfc)) = table2array(Diet(tfr,tfc));
A.dc(isnan(A.dc)) = 0;

% A.dc(1:nr,1:nc-1) = table2array(Diet(1:nr,1:nc-1));
% A.import   = 1 - sum(A.dc)'; %Diet.Outside;

% Flow to detritus

% cname = Base.Properties.VariableNames;
% tmp = regexp(Base.Properties.VariableDescriptions, 'Original column heading: ''(.*)''', 'tokens', 'once');
% isemp = cellfun('isempty', tmp);
% tmp = cat(1, tmp{:});
% cname(~isemp) = tmp;

ndet = A.ngroup - A.nlive;

idxdet = (1:ndet)+10;
idxfish = (1:A.ngear)+max(idxdet);
idxdisc = (1:A.ngear)+max(idxfish);

% detname = Base.Group(isdet);
% [~, detloc] = ismember(detname, bcol);
A.df = table2array(Base(isgroup,idxdet));

% Fisheries

A.fleet = Base.Group(isgear);
% [~, fleetloc] = ismember(A.fleet, bcol);
% [tf, discloc] = ismember(cellfun(@(x) ['disc_' x], A.fleet, 'uni', 0), bcol);
% if ~any(tf)
%     [tf, discloc] = ismember(cellfun(@(x) [x '_disc'], A.fleet, 'uni', 0), bcol);
% end
% if ~any(tf)
%     error('Discard columns not found');
% end

A.landing = table2array(Base(isgroup,idxfish));
A.discard = table2array(Base(isgroup,idxdisc));

A.discardFate = table2array(Base(isgear,idxdet));

% Stanza data

[jtf, jloc] = ismember((1:A.ngroup)', Juvs.JuvNum);
[atf, aloc] = ismember((1:A.ngroup)', Juvs.AduNum);
    
A.stanza = zeros(A.ngroup,1);
A.stanza(jtf) = jloc(jtf);
A.stanza(atf) = aloc(atf);

A.ageStart = zeros(A.ngroup,1);
A.ageStart(atf) = Juvs.RecAge*12;

A.vbK = ones(A.ngroup,1) * -1;
A.vbK(jtf) = Juvs.VonBK;
A.vbK(atf) = Juvs.VonBK;

bab = A.ba./A.b;
babtmp = [bab(Juvs.JuvNum) bab(Juvs.AduNum)];
if all(babtmp(:,1) == babtmp(:,2))
    bab = babtmp(:,1);
else
    bab = num2cell(babtmp,2);
end

nstanza = size(Juvs,1);
A.stanzadata = dataset(...
    {(1:size(Juvs,1))', 'StanzaID'}, ...
    {Juvs.StanzaName,   'StanzaName'}, ...
    {nan(nstanza,1),    'HatchCode'}, ...
    {bab,               'BABsplit'}, ... 
    {nan(nstanza,1),    'WmatWinf'}, ...
    {Juvs.RecPower,     'RecPower'}, ...
    {nan(nstanza,1),    'FixedFecundity'}, ...
    {nan(nstanza,1),    'LeadingLifeStage'}, ...
    {nan(nstanza,1),    'EggsAtSpawn'}, ...
    {nan(nstanza,1),    'LeadingCB'});  % Note: mostly placeholders right now


% Pedigree (At the moment, Rpath allows pedigree values for B, PB, QB, DC,
% and each gear's catch, while I focus on B, PB, QB, DC, EE, and GE)

A.pedigree = nan(A.ngroup, 6);
A.pedigree(:,1) = Ped.B;
A.pedigree(:,2) = Ped.PB;
A.pedigree(:,3) = Ped.QB;
A.pedigree(:,4) = Ped.Diet;

function name = parsecolname(T)
name = T.Properties.VariableNames;
tmp = regexp(T.Properties.VariableDescriptions, 'Original column heading: ''(.*)''', 'tokens', 'once');
isemp = cellfun('isempty', tmp);
tmp = cat(1, tmp{:});
name(~isemp) = tmp;
