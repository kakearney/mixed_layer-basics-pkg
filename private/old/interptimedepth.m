function [newdata, str] = interptimedepth(data, Grd)
%INPTERPTIMEDEPTH Interpolate time vs depth input data

dv = data(2:end, 1:6);
z = data(1,7:end);
olddata = data(2:end,7:end);

% Replicate over years if climatology is given

if length(unique(dv(:,1))) == 1
    yrs = Grd.start_date(1):Grd.end_date(1);
    nyr = length(yrs);
    nperyear = size(olddata,1);
    dvtemp = repmat(dv, nyr, 1);
    dvtemp(:,1) = kron(yrs', ones(nperyear,1));
    dv = dvtemp;
    olddata = repmat(olddata, nyr, 1);
end

% Interpolate

t = (datenum(dv) - datenum(Grd.start_date))*86400;

str = findgap(Grd.time, t, Grd.start_date);
gapsfound = ~isempty(str);

newdata = interp2(t, z, olddata', Grd.time, Grd.z);

if gapsfound
    
    missingt = all(isnan(newdata), 1);
    
    if all(missingt)
        data1 = interp1(z, olddata(1,:), Grd.z);
        newdata = repmat(data1, 1, size(newdata,2));
%         error('Can''t have sim finish before data starts, or vice versa');
    else
  
        first = find(~missingt, 1, 'first');
        last = find(~missingt, 1, 'last');
        [nrow,ncol] = size(newdata);
        before = (1:ncol) < first;
        after = (1:ncol) > last;


        nbefore = sum(before);
        nafter = sum(after);

        newdata(:, missingt & before) = newdata(:,first) * ones(1,nbefore);
        newdata(:, missingt & after)  = newdata(:,last)  * ones(1, nafter);

        missingz = all(isnan(newdata), 2);
        first = find(~missingz, 1, 'first');
        last = find(~missingz, 1, 'last');
        [nrow,ncol] = size(newdata);
        before = (1:nrow)' < first;
        after = (1:nrow)' > last;

        nbefore = sum(before);
        nafter = sum(after);

        newdata(missingz & before,:) = ones(nbefore,1) * newdata(first,:);
        newdata(missingz & after,:)  = ones(nafter,1)  * newdata(last,:);
    end
    
end