function Biovars = setdiapauseparams(BioIn, Biovars, Grd)

% Params that need to be tranferred from ZL to ZL1 and Zl2

if BioIn.isnem
    epvars = {'gs', 'ge', 'mor0', 'Kmor'};
else
    epvars = {'gs', 'ge', 'm0exp'};
end
for ii = 1:length(epvars)
    Biovars.(epvars{ii})([Biovars.idx.zl1 Biovars.idx.zl2]) = Biovars.(epvars{ii})(Biovars.idx.zl);
end

% Mark each time step as either a swim up (1), swim down (-1), or
% no-directed-movement (0) one for the diapause group

yr = Grd.start_date(1):Grd.end_date(1);

d1 = datevec(BioIn.dpEnd);
d2 = datevec(datenum(d1) + BioIn.dpEndSpan);
d3 = datevec(BioIn.dpStart);

bin = [];
val = [];
for iy = 1:length(yr)
    bintmp = [d1; d2; d3];
    bintmp(:,1) = yr(iy);
    bin = [bin; bintmp];
    val = [val; 1; 0; -1];
end
if datenum(Grd.start_date) < datenum(bin(1,:))
    bin = [Grd.start_date; bin];
    val = [0; val];
end
if datenum(Grd.end_date) > datenum(bin(end,:));
    bin = [bin; Grd.end_date];
    val = [val; val(end)];
end
bin = datenum(bin);

dnsim = datenum(Grd.start_date) + Grd.time/86400;
[n, binidx] = histc(dnsim, bin);

Biovars.zlswim = val(binidx);
idx = find(Biovars.zlswim < 0, 1);
Biovars.zlswim(1:idx-1) = 0;

% Set when they swim up and down

% dv = datevec(datenum(Grd.start_date) + Grd.time/86400);
% day = datenum(dv) - datenum([dv(:,1) ones(Grd.nt,2)]);
% 
% Biovars.zlswim = zeros(Grd.nt,1);
% Biovars.zlswim(day >= BioIn.dday(1) & day <= BioIn.dday(2)) = 1; % swim up
% Biovars.zlswim(day >= BioIn.dday(3) | day <  BioIn.dday(1)) = -1; % stay down

% Mark which time steps involve transfer between ZL2 and ZL1.  zlsplit hols
% the fraction of the ZL1 population to be transfered to ZL2 at each time
% step.  zlcombine is a logical marker indicating when to recombine the
% two.

d1 = datevec(BioIn.dpStart);
d2 = datevec(datenum(d1) + BioIn.dpSpan);

dt = mean(diff(dnsim)); % time step in days
nstep = BioIn.dpSpan./dt;

frac = BioIn.dpPercent/100;

bin = [];
val = [];
for iy = 1:length(yr)
    bintmp = [d1; d2];
    bintmp(:,1) = yr(iy);
    bin = [bin; bintmp];
    val = [val; frac; 0];
end
if datenum(Grd.start_date) < datenum(bin(1,:))
    bin = [Grd.start_date; bin];
    val = [0; val];
end
if datenum(Grd.end_date) > datenum(bin(end,:));
    bin = [bin; Grd.end_date];
    val = [val; val(end)];
end
bin = datenum(bin);

frac = BioIn.dpPercent; 
dnsim = datenum(Grd.start_date) + Grd.time/86400;
[n, binidx] = histc(dnsim, bin);

Biovars.zlsplit = val(binidx);
Biovars.zlcombine = [false; Biovars.zlswim(2:end) == 0 & Biovars.zlswim(1:end-1) == 1];

% % Set when to split and recombine the two groups
% 
% Biovars.zlsplit   = [false; Biovars.zlswim(2:end) == -1 & Biovars.zlswim(1:end-1) ~= -1];
% d1 = find(Biovars.zlsplit, 1);
% Biovars.zlswim(1:d1-1) = 0; % keep as one until first down day
% Biovars.zlcombine = [false; Biovars.zlswim(2:end) == 0 & Biovars.zlswim(1:end-1) == 1];

% Some extras

Biovars.t = Grd.time;
% Biovars.dfrac = BioIn.dfrac;

