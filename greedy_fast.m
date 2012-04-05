function [X,done] = greedy_fast( X , PROB )

% if mod(X.N_cones , 100) == 99
%     fprintf('   LL min %f max %f',min(X.greedy_ll(isfinite(X.greedy_ll))), max(X.greedy_ll(:)))
%     figure(1)
%     pe = plotable_evidence(X.greedy_ll) ; 
%     imagesc(pe)
% 
%     figure(2)
%     plot_cones(X.state,PROB) ;
%     'word'
% end

p = PROB.LL-1e20*repmat(X.excluded,[1 1 3]) ;
[mm,I] = max(p(:)) ;
[mx,my,mc] = ind2sub(size(PROB.LL),I) ;

if mm>0
    done = false ;
        
    newX = flip_LL( X , [mx my mc] , PROB , [1 1]) ;
    if newX.ll>=X.ll
        X = update_X({newX},1,false) ;
    else
        [~,indices] = not_excluded( X, mx, my ) ;
        X.excluded(indices) = true ;
    end
    
    try
        fprintf('   #keep_cones %d, #keep_GCs %d    mm-dll %f   mm-PROB.ll %f',...
            nnz(X.keep_cones),numel(X.keep_GCs),mm-newX.ll+oldll,mm-PROB.LL(mx,my,mc)) ;
    end
else
    done = true ;
end

end