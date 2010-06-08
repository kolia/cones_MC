function X = update_X(oldX,X)

check_X(oldX)

deleted = X.diff.deleted(:,1) + (X.diff.deleted(:,2)-1)*X.M0 ;

for d=1:4
    X.contacts{d}(deleted,:) = false ;
    X.contacts{d}(:,deleted) = false ;
end

for d=1:4
    for i=1:size(X.diff.added,1)
        X = make_contacts(X , X.diff.added(i,1) , X.diff.added(i,2)) ;
    end
end

check_X(X)

end