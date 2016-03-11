function mort = nonpredmort(nemflag, bv, A, fe2n, nb, nz)

mort = zeros(nb+2,nb+2,nz);

% N

if nemflag

    [mor0, temp] = meshgrid(A.mor0, A.temp);
    kmor = repmat(A.Kmor', size(mor0,1), 1);
    mortn = bsxfun(@times, bv.^2, tempdep(mor0, kmor, temp));
    mort(1:nb,A.idx.pon,:) = permute(mortn, [2 3 1]);

else
    
    mortn = bsxfun(@times, bsxfun(@power, bv, A.m0exp'), A.m0coef');
    mort(1:nb,A.idx.pon,:) = mortn'; % N flux (per critter -> PON)

end

% Si

mort(A.idx.plsi,A.idx.opal,:) = mortn(:,A.idx.pl) .* A.RSiN;

% Fe

mortfe = mortn(:,[A.idx.ps A.idx.pl]) .* fe2n(:, [A.idx.ps A.idx.pl]);
mort([A.idx.psfe A.idx.plfe],A.idx.pofe,:) = permute(mortfe, [2 3 1]);

%------------------------
% Q10-style temperature 
% dependence
%------------------------

function td = tempdep(a, b, temp)
td = a .* exp(b.*temp);