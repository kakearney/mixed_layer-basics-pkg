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

% Set when they swim up and down

dv = datevec(datenum(Grd.start_date) + Grd.time/86400);
day = datenum(dv) - datenum([dv(:,1) ones(Grd.nt,2)]);

Biovars.zlswim = zeros(Grd.nt,1);
Biovars.zlswim(day >= BioIn.dday(1) & day <= BioIn.dday(2)) = 1; % swim up
Biovars.zlswim(day >= BioIn.dday(3) | day <  BioIn.dday(1)) = -1; % stay down

% Set when to split and recombine the two groups

Biovars.zlsplit   = [false; Biovars.zlswim(2:end) == -1 & Biovars.zlswim(1:end-1) ~= -1];
d1 = find(Biovars.zlsplit, 1);
Biovars.zlswim(1:d1-1) = 0; % keep as one until first down day
Biovars.zlcombine = [false; Biovars.zlswim(2:end) == 0 & Biovars.zlswim(1:end-1) == 1];

% Some extras

Biovars.t = Grd.time;
Biovars.dfrac = BioIn.dfrac;

