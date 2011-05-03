function replay_select_X(cone_map,dX,T,which_ones,do_something)
% cone_map can be the initial cone_map produced by exact_LL_setup.m (or
% equivalently , by main_CAST.m), or it can be the final cone_map after 
% a run, either way.
% dX can be found as a field of X.
% which_ones is a list of iteration numbers which we would like to
% do_something with.  do_something will be called with arguments X,
% cone_map, and the iteration number i.
% You can have do_something save or plot X at which_ones iterations.

f = @do_it ;
function i = do_it(i,X)
    if ismember(i,which_ones)
        do_something(X,cone_map,i) ;
    end
    i = i+1 ;
end

dX = dX(1:max(which_ones)+5,:) ;

fold_X( cone_map.initX , dX , cone_map , T , 0 , f ) ;

end