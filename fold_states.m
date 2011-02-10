function r = fold_states( state , dX , r , f )
% r = fold_states( state , dX , r , f )
% 
% Replay MCMC sequence by starting with initial state 'state' and updating
% it using dX.  After each update, the function f is applied to r and
% state.  fold_states folds f over the MCMC state sequence defined by
% (state,dX), starting with initial value r.  Also returns final state.


samples = size(dX,1)  ;
[N,M]   = size(state) ;

for i=1:samples

    % Cone deletions
    inds = find( (dX(i,3:3:end) == 0) & (dX(i,1:3:end) ~= 0) ) ;
    x = full( dX(i,3*(inds-1)+1) ) ;
    y = full( dX(i,3*(inds-1)+2) ) ;
    state(x+N*(y-1)) = 0 ;

    % Cone additions
    for c=1:3
        inds = find( dX(i,3:3:end) == c ) ;
        x = dX(i,3*(inds-1)+1) ;
        y = dX(i,3*(inds-1)+2) ;
        state(x+N*(y-1)) = c ;
    end
    
    r = f(r,state) ;
    
end

end