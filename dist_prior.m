function dll = dist_prior(X,i,j,type,e,sigma,N_cones_prior)

if strcmp( type , 'color_change' )
    dll = 0 ;
    return
end

[M0,M1] = size(X.state) ;

D = ceil(e) ;

sx      = max(1,i-D):min(M0,i+D) ;
sy      = max(1,j-D):min(M1,j+D) ;

[rx,ry] = find( X.state(sx,sy) > 0 ) ;
rx = sx(1) + rx - 1 ;
ry = sy(1) + ry - 1 ;

dll = sum( sigma( (rx-i).^2 + (ry-j).^2 ) ) + N_cones_prior(X.N_cones) ;

if strcmp( type , 'deletion' )
    dll = - dll ;
end

end