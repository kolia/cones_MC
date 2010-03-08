function swap = swap_LL( swap , flips , LL , withLL )
% Simply applies 'flips' to both systems in a swapper, and exposes swap.ll
% as the resulting joint log-likelihood and swap.state, the state of swap.X
% to allow for accumulation of observables by flip_MCMC.

swap.X      =     LL( swap.X    , flips ) ;
swap.with   = withLL( swap.with , flips ) ;

swap.ll = (swap.X.ll + swap.with.ll) ;
swap.state = swap.X.state ;

end