%% Eecho: Volcano again, now with new concentration-dependent iron 
%  scavenging, realistic interannual forcing, and limited pelagic predation 
           
% Create food webs

rerun = false;

if rerun

    Tmp = load('eraInterimPhysics');
    Era  = Tmp.Era;  % full timeseries
    
    EraS = Tmp.EraS; % sample-year
    EraS.syear = 1900; % b/c foodwebsetup expects this
    EraS.eyear = 1902; 

    Bgc = esabgcsetup;

    % Diapause on, with no predation at depth

    Bgc.diapause = true;
    Bgc.dfrac = 0.9;
%     Bgc.predatdepth = false;

    % Single-food web arrays, used to compare the more- and less-simplified
    % versions of the food web.

    FwSingle(1) = esafoodwebsetupfull(EraS, Bgc, 1, ...
             'runnemuro', false, ...
             'nemfile', 'eechoNemforcal.nc', ...
             'cvtparams', {'wwCfrac', 0.03});
    FwSingle(2) = esafoodwebsetup(EraS, Bgc, 1, ...
             'runnemuro', false, ...
             'nemfile', 'eechoNemforcal.nc', ...
             'cvtparams', {'wwCfrac', 0.03});
    
    Fw = esafoodwebsetupfull(EraS, Bgc, 100, ...
             'runnemuro', false, ...
             'nemfile', 'eechoNemforcal.nc', ...
             'cvtparams', {'wwCfrac', 0.03});

    save eechosetup EraS Era Bgc Fw FwSingle;
else
    Tmp = load('eechosetup');
    EraS = Tmp.EraS;
    Era = Tmp.Era;
    Bgc = Tmp.Bgc;
    Fw  = Tmp.Fw;
    FwSingle = Tmp.FwSingle;
end

%-----------------
% Spinup/baseline 
% sim
%-----------------

% Adjust iron scavenging (not relevant for baseline iron levels, but might
% as well set this up ahead of time)

Bgc.scavthresh = 0.6e-6; % Johnson et al. 1997 theshold
Bgc.scavfac = 0.08;      % Tuned to get 3-5 day halflife

In = formatforwce(Era, Bgc, Fw, 'constant');   

% Most bird and mammal non-pred losses leave the system (i.e. to land or 
% deep water), so reroute these losses from PON to out.  Also, prescribe
% 90% of non-predatory losses as fisheries losses for commercially-fished
% critters.  

[sname, lname] = readwcevarnames('/Volumes/MyPassportKak/wceSims/ccharlie/ccharlie_base/ccharlie_base001.nc');

bigguys = {...
    'Albatross'
    'Sperm whales'
    'Toothed whales'
    'Elephant seals'
    'Seals,dolphins'
    'Fulmars'
    'Skuas,Jaegers'
    'Puffins,Shearwaters,Storm Petrels'
    'Kittiwakes'
    'Sharks'};

fishedguys = {...
    'Pomfret'
    'Chum salmon'
    'Chinook,coho,steelhead'
    'Sockeye,Pink'
    'Saury'
    'Sergestid shrimp'};
    
[tf1,loc1] = ismember(bigguys, lname);
[tf2,loc2] = ismember(fishedguys, lname);

n1 = length(loc1);
n2 = length(loc2);

newreroute = cell(n1+n2,5);
[newreroute{:,1}] = deal('mor');
newreroute(:,2) = sname([loc1; loc2]);
[newreroute{:,3}] = deal('PON');
[newreroute{1:n1,4}] = deal('out');
[newreroute{n1+1:end,4}] = deal('fish');
[newreroute{:,5}] = deal(0.9);
newreroute = [In(1).reroute; newreroute];

% Limit pelagic fishes and seabirds to only feed on the upper water column,
% allowing copepods to escape them when they enter diapause.  All other
% critters remain able to prey on the entire water column.

zvis = [0 250 400 1000];
pelvis = [1 1 0 0];

pelagics = {...
    'Albatross'
    'Fulmars'
    'Chinook,coho,steelhead'
    'Skuas,Jaegers'
    'Pomfret'
    'Puffins,Shearwaters,Storm Petrels'
    'Kittiwakes'
    'Sockeye,Pink'
    'Pelagic forage fish'
    'Saury'
    'Chum salmon'};

ispelagic = ismember(Fw.EweinEns(1).name, pelagics);

vis = ones(length(zvis),Fw.EweinEns(1).ngroup);
vis(:,ispelagic) = pelvis' * ones(1,sum(ispelagic));

preyvis = [-zvis' vis];

% Want propagation to higher trophic levels, so set linear mortality for
% all nekton.  Put zooplankton at 1.5, as compromise between flow up and
% model stability.

m0exp = ones(In(1).Ewein.nlive,1);
m0exp(strncmp(In(1).types(1:In(1).Ewein.nlive), 'z', 1)) = 1.5;
m0exp(strncmp(In(1).types(1:In(1).Ewein.nlive), 'p', 1)) = 2;

% Add these adjustments to all ensemble members

[In(:).m0exp]   = deal(m0exp);
[In(:).reroute] = deal(newreroute);
[In(:).preyvis] = deal(preyvis);
[In(:).verbose] = deal(false); 
     
% Save 2008+ years as baseline comparison.  Also save a start file at 2008
% for use in any perturbation runs.  For testing, save monthly-averages
% throughout.

archparams = {'beginarchive', [NaN datenum(2008,1,1)], ...
              'endarchive',   [NaN NaN], ...
              'tarch',        [-1 86400], ...
              'outputextension', {'monthly', '2008to2010daily'}};

for ii = 1:length(In)
    In(ii).hotstartdn = datenum(2008,1,1);
    In(ii).hotstartsave = fullfile('/Volumes/MyPassportKak/wceSims/eecho/hs', sprintf('hs%03d.mat', ii));
end
if ~exist('/Volumes/MyPassportKak/wceSims/eecho/hs', 'dir')
    mkdir('/Volumes/MyPassportKak/wceSims/eecho/hs');
end

runsim = @(name, A) runmixedlayer(A, ...
    'folder', '/Volumes/MyPassportKak/wceSims/eecho', ...
    'name' , ['eecho_' name], ...
    'idx', 1:100, ...
    'usepar', true, ...
    'postprocess', true, ...
    'consolidatebio', false, ...
    archparams{:});

% notifier('3057244998@vtext.com', @() runsim('base', In));

%-----------------
% Volcano
%-----------------

% Surface iron timeseries: Use constant background value, except during the
% Aug. 9-11, 2008 (approximately when ash cloud crossed ESA region) at
% which point add a spike of 1000X  

nyr = Era.eyear - Era.syear + 1;
tfe = datelist([Era.syear 1 1 0 0 0], [Era.eyear 12 31 0 0 0], 1:12);

feflux = Bgc.feflx * ones(nyr*12,1);

tspike = [...
    2008 8 08 23 0 0 % Make onset as steep as possible
    2008 8 09 0  0 0
    2008 8 11 0  0 0
    2008 8 12 0  0 0];
tspike = datenum(tspike);
fac = ones(size(tspike));
fac(2:end-1) = 1000;
fespike = interp1(tfe, feflux, tspike) .* fac;

Bgc.feflxTs = sortrows([[tfe; tspike], [feflux; fespike]]);

InV = formatforwce(Era, Bgc, Fw, 'timeseries');

% Add the extras

[InV(:).m0exp]   = deal(m0exp);
[InV(:).reroute] = deal(newreroute);
[InV(:).preyvis] = deal(preyvis);
[InV(:).verbose] = deal(false);

% Hotstart from baseline

for ii = 1:length(InV)
    InV(ii).hotstartload = fullfile('/Volumes/MyPassportKak/wceSims/eecho/hs', sprintf('hs%03d.mat', ii));
end

% Check for files whose baseline runs didn't make it to the hotstart time

ismissing = cellfun(@(x) ~exist(x, 'file'), {InV.hotstartload});

% Run (with a slightly different caller, since need to supply
% expected-missing runs for postprocessing)

runsimhs = @(name, A) runmixedlayer(A, ...
    'folder', '/Volumes/MyPassportKak/wceSims/eecho', ...
    'name' , ['eecho_' name], ...
    'idx', 1:100, ...
    'usepar', true, ...
    'postprocess', true, ...
    'consolidatebio', false, ...
    'idxmissing', find(ismissing), ...
    archparams{:});

% notifier('3057244998@vtext.com', @() runsimhs('volc', InV));
         
%-----------------
% Volcano
% spring eruption
%-----------------

nyr = Era.eyear - Era.syear + 1;
tfe = datelist([Era.syear 1 1 0 0 0], [Era.eyear 12 31 0 0 0], 1:12);

feflux = Bgc.feflx * ones(nyr*12,1);

tspike = [...
    2008 3 08 23 0 0 % Make onset as steep as possible
    2008 3 09 0  0 0
    2008 3 11 0  0 0
    2008 3 12 0  0 0];
tspike = datenum(tspike);
fac = ones(size(tspike));
fac(2:end-1) = 1000;
fespike = interp1(tfe, feflux, tspike) .* fac;

Bgc.feflxTs = sortrows([[tfe; tspike], [feflux; fespike]]);

InVearly = formatforwce(Era, Bgc, Fw, 'timeseries');

% Add the extras

[InVearly(:).m0exp]   = deal(m0exp);
[InVearly(:).reroute] = deal(newreroute);
[InVearly(:).preyvis] = deal(preyvis);
[InVearly(:).verbose] = deal(false);

% Hotstart from baseline

for ii = 1:length(InV)
    InVearly(ii).hotstartload = fullfile('/Volumes/MyPassportKak/wceSims/eecho/hs', sprintf('hs%03d.mat', ii));
end

% notifier('3057244998@vtext.com', @() runsimhs('early', InVearly));

%-----------------
% Volcano
% summer eruption
%-----------------

nyr = Era.eyear - Era.syear + 1;
tfe = datelist([Era.syear 1 1 0 0 0], [Era.eyear 12 31 0 0 0], 1:12);

feflux = Bgc.feflx * ones(nyr*12,1);

tspike = [...
    2008 6 08 23 0 0 % Make onset as steep as possible
    2008 6 09 0  0 0
    2008 6 11 0  0 0
    2008 6 12 0  0 0];
tspike = datenum(tspike);
fac = ones(size(tspike));
fac(2:end-1) = 1000;
fespike = interp1(tfe, feflux, tspike) .* fac;

Bgc.feflxTs = sortrows([[tfe; tspike], [feflux; fespike]]);

InVsummer = formatforwce(Era, Bgc, Fw, 'timeseries');

% Add the extras

[InVsummer(:).m0exp]   = deal(m0exp);
[InVsummer(:).reroute] = deal(newreroute);
[InVsummer(:).preyvis] = deal(preyvis);
[InVsummer(:).verbose] = deal(false);

% Hotstart from baseline

for ii = 1:length(InV)
    InVsummer(ii).hotstartload = fullfile('/Volumes/MyPassportKak/wceSims/eecho/hs', sprintf('hs%03d.mat', ii));
end

% notifier('3057244998@vtext.com', @() runsimhs('summer', InVsummer));

%-----------------
% Volcano 
% no diapause
%-----------------

InVnoD = InV;
[InVnoD.dfrac] = deal(0); % 

notifier('3057244998@vtext.com', @() runsimhs('nodiapause', InVnoD));

%-----------------
% Volcano
% delayed diapause
%-----------------

InVDelayDiapause = InV;
[InVDelayDiapause.dday] = deal(InV(1).dday + [0 0 21]);

notifier('3057244998@vtext.com', @() runsimhs('delaydiapause', InVDelayDiapause));

%-----------------
% Perpetually-
% fertilized
%-----------------

% Tried using the same hot-start, but this leads to too much drift in the
% big guys, who need to adjust slowly to the new iron levels.  Need to run
% with a full spinup so the drift settles out by 2008. 

BgcHi = Bgc;
BgcHi.feflx = Bgc.feflx * 1000;
InHi = formatforwce(Era, BgcHi, Fw, 'constant');

[InHi(:).m0exp]   = deal(m0exp);
[InHi(:).reroute] = deal(newreroute);
[InHi(:).preyvis] = deal(preyvis);
[InHi(:).verbose] = deal(false);

notifier('3057244998@vtext.com', @() runsim('constanthi', InHi));

%-----------------
% Volcano
% no predation at
% depth
%-----------------

vis = pelvis' * ones(1,Fw.EweinEns(1).ngroup);
preyvisNone = [-zvis' vis];

InVnoPred = InV;
[InVnoPred(:).preyvis] = deal(preyvisNone);

% notifier('3057244998@vtext.com', @() runsimhs('nopred', InVnoPred));

%-----------------
% Volcano
% all prey at 
% depth
%-----------------

vis = ones(length(zvis),Fw.EweinEns(1).ngroup);
preyvisAll = [-zvis' vis];

InVallPred = InV;
[InVallPred(:).preyvis] = deal(preyvisAll);

% notifier('3057244998@vtext.com', @() runsimhs('allpred', InVallPred));