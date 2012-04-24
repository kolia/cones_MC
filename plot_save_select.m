function plot_save_select(cone_map,dX,which_ones)
% cone_map can be the initial cone_map produced by exact_LL_setup.m (or
% equivalently , by main_CAST.m), or it can be the final cone_map after 
% a run, either way.
% dX can be found as a field of X.
% which_ones is a list of iteration numbers which we would like to
% do_something with.

T = [1 1] ;

do = @plot_save ;
    function plot_save(X,PROB,i)
        save(sprintf('X_iteration_%d',i),'X')
        plot_cones_matlab({X},PROB) ;
        title(sprintf('Iteration %d   LL%8d',i,floor(X.ll)),'FontSize',22)
        drawnow
        saveas(gcf,sprintf('X_at_iteration_%d',i),'fig')        
    end

replay_select_X(cone_map,dX,T,which_ones,do)

end

