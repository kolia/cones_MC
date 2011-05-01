function X = greedy( X , PROB , update_X )

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

old_ll  = X.ll ;
ll      = old_ll * ones(M0,M1,PROB.N_colors) ;

for x=1:M0
    for y=1:M1
        % propose addition of new cone of each color
        for c=1:PROB.N_colors
            if X.state(x,y) ~= c
                sample = change_cone( X , [x y c] , PROB , [1 1]) ;
                ll(x,y,c) = sample.ll ;
            end
        end
    end
end

m = max(ll(:)) ;

if m>old_ll
    [mx,my] = find(ll == m) ;
    
    mx = mx(1) ;
    my = my(1) ;
    
    mc = 1+floor((my-1)/M1) ;
    my = 1+mod(my-1,M1) ;
    
    X = change_cone( X , [mx my mc] , PROB , [1 1]) ;
    X = update_X({X}) ;
end

end