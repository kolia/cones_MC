function ll = calculate_LL( X , PROB , T )

if X.N_cones>0
    ll  = 0.5 * T(1) * (sum( PROB.N_cones_terms .* sum(abs(X.sparse_STA_W_state)>0,2)) + ...
            sum( (PROB.quad_factors .* X.contributions).^T(2) ) ) ;
else
    ll = 0 ;
end

end