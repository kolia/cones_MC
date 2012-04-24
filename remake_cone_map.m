function cone_map = remake_cone_map( X )
% Use information in X to reconstruct cone_map.

load(sprintf('%s/stas', X.datafolder))
load(sprintf('%s/cone_params', X.datafolder))
cone_map = transfer_info( X, struct ) ;
cone_map = exact_LL_setup(stas,cone_params,cone_map) ;

end