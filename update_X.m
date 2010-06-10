function X = update_X(oldX,X)

check_X(oldX)

if ~isempty(X.diff.deleted)
    deleted = X.diff.deleted(:,1) + (X.diff.deleted(:,2)-1)*X.M0 ;
    
    for d=1:4
        X.contact{d}(deleted,:) = false ;
        X.contact{d}(:,deleted) = false ;
    end
end


for d=1:4
    for i=1:size(X.diff.added,1)
        X = make_contacts(X , X.diff.added(i,1) , X.diff.added(i,2)) ;
    end
end

check_X(X)

X.diff.added = [] ;
X.diff.deleted = [] ;

end