function make_plots(greed,mcmc,cast,cone_map)

if nargin<4
    cone_map = remake_cone_map(greed.initX) ;
end

folder_name = cone_map_string(cone_map);

filename = ['plots_for_' folder_name] ;
mkdir( filename )

here = pwd ;

cd(filename)

[~,ll] = get_best( cast ) ;

try load(['../confident_' folder_name])
catch
    selector = @(n) (n>10000) && (mod(n,20) == 0) ;
    dX = cast{find(ll==max(ll),1)}.X.dX ;
    confident = confident_cones( greed.X , dX , cone_map , selector ) ;
end
save(['../confident_' folder_name], 'confident') ;
plot_cone_field( confident , cone_map )

plot_LL_ncones( greed , mcmc , cast , cone_map )
plot_Greedy_MCMC_CAST( greed , mcmc , cast , cone_map.NICE )


% [sta,invww] = denoised_sta( greed.initX , cast{1}.X.dX , cone_map, selector ) ;
% make_sta_plots( sta , invww, 'denoised' )

cd(here)

end