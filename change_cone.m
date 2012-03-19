function X = change_cone( X , a , PROB , T )

% check_X(X)

for i=1:size(a,1)
    x = a(i,1) ;
    y = a(i,2) ;
    c = a(i,3) ;
    free = ~c || ~X.state(x,y) ;
    
    if  c && ~X.state(x,y)  % check exclusion before adding cone
        free = not_excluded(X,x,y) ;
    end

    if free                 % add or delete cone
        X    = flip_LL( X , [x y c] , PROB , T ) ;
    else
        X.ll = -Inf ;
        break
    end
    
end

end