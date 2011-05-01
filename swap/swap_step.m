function [X1,X2] = swap_step(X1,T1,X2,T2,cone_map)

swapX = swap_closure( X1,T1, X2,T2 , cone_map) ;
swapX = flip_MCMC( swapX{1}, swapX(2:end), @update_swap ) ;
X1  = swapX.X ;
X2  = swapX.with ;

end
