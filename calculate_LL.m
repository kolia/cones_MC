function ll = calculate_LL( X , PROB , T )

if X.N_cones>0
    ll  = 0.5 * T(1) * (sum( PROB.N_cones_terms .* sum(X.sparse_STA_W_state>0,2)) + ...
            sum(X.contributions)) ;
else
    ll = 0 ;
end

end