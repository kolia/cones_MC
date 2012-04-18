function [X,done] = greedy_hot( X , PROB )

% p = log-posterior increase by adding cones at all locations/colors
p = PROB.LL-1e20*repmat(X.excluded,[1 1 3]) ;   % X.excluded = hard prior

% best cone addition location and color
[LL_increase,I] = max(p(:)) ;
[mx,my,mc] = ind2sub(size(PROB.LL),I) ;

if LL_increase>0
    done = false ;

    % calculate true LL of new configuration (with no approximation)
    newX = flip_LL( X , [mx my mc] , PROB , [1 1]) ;

    % if the LL did in fact increase, add cone
    if newX.ll>=X.ll
        X = update_X({newX},1,PROB.track_contacts) ;
    % otherwise, poor approximation: exclude region around proposed cone
    else
        [~,indices] = not_excluded( X, mx, my ) ;
        X.excluded(indices) = true ;
    end
    
    % DISPLAY to stdout
    try
        fprintf('   #keep_cones %d, #keep_GCs %d    mm-dll %f   mm-PROB.ll %f',...
            nnz(X.keep_cones),numel(X.keep_GCs),mm-newX.ll+oldll,mm-PROB.LL(mx,my,mc)) ;
    end
else
    done = true ;
end

end