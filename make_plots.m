function make_plots(greed,mcmc,cast,cone_map)

if nargin<4
    cone_map = remake_cone_map(greed.initX) ;
end

filename = ['plots_for_' cone_map_string(cone_map)] ;
mkdir( filename )

here = pwd ;

cd(filename)

plot_LL_ncones( greed , mcmc , cast , cone_map )
plot_Greedy_MCMC_CAST( greed , mcmc , cast , cone_map.NICE )

% [sta,invww] = denoised_sta( greed.initX , cast{1}.X.dX , cone_map ) ;
% 
% make_sta_plots( sta , invww, 'denoised' )

cd(here)

end