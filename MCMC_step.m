function X = MCMC_step(X,PROB,T)
LL = @(x)get_LL(x,PROB,T) ;
X  = flip_MCMC( X, move( X, PROB , T), @update_X, LL ) ;
end
