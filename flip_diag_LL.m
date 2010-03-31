function X = flip_diag_LL( X , x , y , c , LL , N_cones_factor , dprior )


[M0,M1,N_colors] = size(LL) ;

% initialize
if x < 1 || y < 1
    X.data_ll   = 0 ;
    X.prior     = 0 ;
    X.ll        = 0 ;
    X.N_cones   = 0 ;

    X.state     = sparse(zeros(M0,M1)) ;
%     X.state     = zeros(M0,M1) ;
    
% update X with flip
else    
    old_color     = X.state(x,y) ;
    
    if old_color
        if c == old_color || c == 0
            type = 'deletion' ;
            X.data_ll = X.data_ll - LL(x,y,old_color) + N_cones_factor ;
            X.state(x,y) = 0 ;
            X.N_cones = X.N_cones - 1 ;
        else
            type = 'color_change' ;
            X.data_ll = X.data_ll + LL(x,y,c) - LL(x,y,old_color) ;
            X.state(x,y) = c ;
        end
    else
        type = 'addition' ;
        X.data_ll = X.data_ll + LL(x,y,c) - N_cones_factor ;
        X.state(x,y) = c ;
        X.N_cones = X.N_cones + 1 ;
    end
    
    X.prior  = X.prior + dprior(X,x,y,type) ;
    
    X.ll      = X.data_ll + X.prior ;
        
end