function X = initialize_X(LL,N_cones_factor,D)

[M0,M1,N_colors] = size(LL) ;

X.masks     = make_masks(D) ;
X.M0        = M0 ;
X.M1        = M1 ;
X.N_colors  = N_colors ;
X.N_cones   = 0  ;
X.N_cones_factor = N_cones_factor ;
X.ll        = 0 ;

% upper bound on anticipated number of cones
X.maxcones  = ceil( M0*M1 * 0.015 ) ;

% sparse int matrix, representing cone positions and colors
X.state     = sparse([],[],[],M0,M1,X.maxcones) ;

% map from cone position to id
X.id        = sparse([],[],[],M0,M1,X.maxcones) ;

% which ids are already assigned to cones
X.taken_ids = false(X.maxcones,1) ;

% sparse int matrix, with number of out-of-border adjacencies
X.outofbounds = sparse([],[],[],M0,M1,2*(M0+M1)) ;
X.outofbounds(:,[1 M1]) = 1 ;
X.outofbounds([1 M0],:) = X.outofbounds([1 M0],:) + 1 ;

% contact forces at four cardinal adjacent positions, indexed by id
for d=1:4
    X.contact{d} = sparse(zeros(X.maxcones)) ;
end

end