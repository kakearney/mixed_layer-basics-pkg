function [newdata, str] = interptime(t, data, Grd)
%INTERPTIME Interpolate input data vs time

str = findgap(Grd.time, t, Grd.start_date);
gapsfound = ~isempty(str);

newdata = interp1(t, data, Grd.time, 'linear');

if gapsfound
    newdataextrap = interp1(t, data, Grd.time, 'nearest', 'extrap');
    isnull = isnan(newdata);
    newdata(isnull) = newdataextrap(isnull);
end