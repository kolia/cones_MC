function ll = calculate_LL( X , PROB , T )

% M0 = PROB.M0 * PROB.SS ;
% M1 = PROB.M1 * PROB.SS ;

if X.N_cones>0
%     invWW = X.invWW ;
%     invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;
%     invWW = sparse(invWW) ;

%     STA_W_state = (X.STA_W_state-PROB.min_STA_W).^T(2) + PROB.min_STA_W ;

    try
        ll = X.ll + 0.5 * T(1) * full( PROB.N_cones_term + ...
            sum( PROB.quad_factor' .* X.dUW_STA.^2)) ;
        
%         ll = 0.5 * T(1) * full( X.N_cones * PROB.N_cones_term + ...
%             sum( PROB.quad_factor' * ((X.WW\STA_W_state')' .* STA_W_state) )) ;
        
%         ll  = 0.5 * T(1) * full( X.N_cones * PROB.N_cones_term + ...
%             sum( PROB.quad_factor' * ((STA_W_state * X.invWW) .* STA_W_state) )) ;
        
%         fprintf('%f,%f,%f,%f\n',X.ll,ll0-X.ll,ll-X.ll,(ll-X.ll)/(ll0-X.ll))
    catch
        'ha'
    end
else
    ll = 0 ;
end

end