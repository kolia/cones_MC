function r = fold_X( X , dX , PROB , T , r , f )
% r = fold_states( state , dX , r , f )
% 
% Replay MCMC sequence by starting with initial state 'state' and updating
% it using dX.  After each update, the function f is applied to r and
% state.  fold_states folds f over the MCMC state sequence defined by
% (state,dX), starting with initial value r.  Also returns final state.

samples = size(dX,1)  ;

for i=1:samples
    inds = find( dX(i,1:3:end) ~= 0 ) ;
    x = full( dX(i,3*(inds-1)+1) ) ;
    y = full( dX(i,3*(inds-1)+2) ) ;
    c = full( dX(i,3*(inds-1)+3) ) ;
    
    X = flip_LL( X , [x(:) y(:) c(:)] , PROB , T ) ;
    r = f(r,X) ;
end

end