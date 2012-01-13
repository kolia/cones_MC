% reap  % load greed, mcmc and cast from hpc cluster

plot_LL_ncones( greed , mcmc , cast )
plot_Greedy_MCMC_CAST( greed , mcmc , cast )

[sta,invww] = denoised_sta( cast{1}.initX , cast{1}.X.dX , cast{1} ) ;

make_sta_plots( sta , invww, 'denoised' )