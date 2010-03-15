function cone_map = greedy_cones( GC_stas , cone_params , cone_map )
%% cone_map = greedy_cones( GC_stas , cone_params )
%  Do greedy search for cone positions.

if nargin<3  ,  cone_map = struct ;   end

% load variables defined in setup_cone_LL
[STA_W,cone_map] = setup_cone_LL(GC_stas , cone_params , cone_map) ;
vars = fields(cone_map) ;
for i=1:length(vars)
    eval(sprintf('%s = cone_map.(vars{i}) ;',vars{i})) ;
end

% define prior
prior_LL = @(X)-1e12*sum(sum( triu(X.overlaps,1) > max_overlap*0.1 )) ;


%% GREEDY SOLUTION
GREED.state  = zeros( 1 , NROI ) ;
best_LL      = -Inf ;

flip_LL      = @(X,flips)flip_color_LL( ...
    X , flips , prior_LL , cell_consts , STA_W' , coneConv , colorDot , sizeROI ) ;


fprintf('\nGREEDY cone finding:\n')
tic
while 1
    GREED = greedy( GREED , flip_LL ) ;
    fprintf('\nCONES:%2d \t \t increase in LL:%f',sum(GREED.state),GREED.ll-best_LL)
    if GREED.ll<=best_LL
        break
    else
        best_LL = GREED.ll ;
    end
    
    GGG = zeros(26,46,3) ;
    GGG(ROI) = GREED.state ;
    imagesc(GGG) ;
    drawnow

end
fprintf('\nGREEDY SOLUTION found in %.1f sec\n',toc)

cone_map.stas        = GC_stas ;
cone_map.cone_RF     = cone_RF ;
cone_map.cone_params = cone_params ;
cone_map.prior_LL    = func2str(prior_LL) ;
cone_map.max_overlap = max_overlap ;
cone_map.GREED       = GREED ;
cone_map.GREED.state = zeros(M0*SS,M1*SS,N_colors) ;
cone_map.GREED.state(ROI) = GREED.state ;

end