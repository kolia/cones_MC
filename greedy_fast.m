function [X,done] = greedy_fast( X , PROB , update_X )

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

[mm,I] = max(PROB.LL(:)) ;
[mx,my,mc] = ind2sub(size(PROB.LL),I) ;

if mm>0
    done = false ;
    
    if newX.ll>=X.ll
        X = update_X({newX},1,false) ;
    end

    try
        fprintf('   #keep_cones %d, #keep_GCs %d    mm-dll %f   mm-PROB.ll %f',...
            nnz(X.keep_cones),numel(X.keep_GCs),mm-newX.ll+oldll,mm-PROB.LL(mx,my,mc)) ;
    end
else
    done = true ;
end

end