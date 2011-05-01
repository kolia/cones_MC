function ll = get_LL( X , PROB , T )
% Get the log-likelihood of X at temperature T.

if T(1) == X.T(1) && T(2) == X.T(2)
    ll = X.ll ;
else    
    ll = calculate_ll( X , PROB , T ) ;
end