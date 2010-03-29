function dll = dist_prior(X,i,j,c,D,sigma,N_cones_prior)

[M0,M1] = size(X.state) ;

sx      = max(1,i-D):min(M0,i+D) ;
sy      = max(1,j-D):min(M1,j+D) ;

[rx,ry] = find( X.state(sx,sy) > 0 ) ;
rx = sx(1) + rx - 1 ;
ry = sy(1) + ry - 1 ;

dists = (rx-i).^2 + (ry-j).^2 ;

dll = - sum( sigma( (rx-i).^2 + (ry-j).^2 ) ) + N_cones_prior(X.N_cones) ;

if ismember(0,dists)
    dll = - dll ;
end

end