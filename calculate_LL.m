function ll = calculate_LL( X , PROB , T )

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

if X.N_cones>0
    [x,y,c] = find(X.state) ;
    X.STA_W_state = PROB.STA_W( x+M0*(y-1)+M0*M1*(c-1) , : )' ; 

    invWW = X.invWW ;
    invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;
    invWW = sparse(invWW) ;

    STA_W_state = (X.STA_W_state-PROB.min_STA_W).^T(2) + PROB.min_STA_W ;

    ll  = 0.5 * T(1) * full( X.N_cones * PROB.N_cones_term + ...
        sum( PROB.quad_factor .* ...
        sum( (STA_W_state * invWW) .* STA_W_state ,2) )) ;
else
    ll = 0 ;
end

end