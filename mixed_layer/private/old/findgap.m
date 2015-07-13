function str = findgap(simlim, datalim, startdate)
%FINDGAP Find time gaps in input data
%
% str = findgap(simlim, datalim, startdate)
%
% Originally part of initialize.m, moved here so could be used by
% biomodules as well.

gap = nan(2,2);
if min(datalim) > min(simlim)
    gap(1,:) = datenum(startdate) + [min(simlim) min(datalim)]./86400;
end
if max(datalim) < max(simlim)
    gap(2,:) = datenum(startdate) + [max(datalim) max(simlim)]./86400;
end
if all(isnan(gap))
    str = '';
elseif isnan(gap(1,1))
    str = [datestr(gap(2,1), 'mm/dd/yyyy') ' - ' datestr(gap(2,2), 'mm/dd/yyyy')];
elseif isnan(gap(2,1))
    str = [datestr(gap(1,1), 'mm/dd/yyyy') ' - ' datestr(gap(1,2), 'mm/dd/yyyy')];
else
    str = [datestr(gap(1,1), 'mm/dd/yyyy') ' - ' datestr(gap(1,2), 'mm/dd/yyyy') ', ' datestr(gap(2,1), 'mm/dd/yyyy') ' - ' datestr(gap(2,2), 'mm/dd/yyyy')];
end