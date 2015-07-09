function bionew = tagalongmixing(bioold, bionew, tagpairs)
%TAGALONGMIXING Mixes biology by keeping ratio between variables constant

for it = 1:size(tagpairs,1)
    
    childold  = bioold(:,tagpairs(it,1));
    parentold = bioold(:,tagpairs(it,2));
    parentnew = bionew(:,tagpairs(it,2));
    
    childnew = parentnew .* childold ./ parentold;
    
    bionew(:,tagpairs(it,1)) = childnew;
    
end
    



