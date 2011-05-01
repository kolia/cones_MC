function [X1,X2] = swap_step(X1,X2,T,PROB)

LL = @(x)get_LL(x.X,PROB,T)+get_LL(x.with,PROB,T) ;
swapX = swap_closure( X1, X2 ,T, PROB) ;
swapX = flip_MCMC( swapX{1}, swapX(2:end), @update_swap , LL ) ;
X1  = swapX.X ;
X2  = swapX.with ;

end
