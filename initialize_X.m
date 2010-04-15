function X = initialize_X(LL,D)

[M0,M1,N_colors] = size(LL) ;

X.masks     = make_masks(D) ;
X.M0        = M0 ;
X.M1        = M1 ;
X.N_colors  = N_colors ;
X.N_cones   = 0  ;
X.N_cones_factor = N_cones_factor ;

% upper bound on anticipated number of cones
X.maxcones  = ceil( M0*M1 * 0.015 ) ;

% sparse int matrix, representing cone positions and colors
X.state     = sparse([],[],[],M0,M1,X.maxcones) ;

% map from cone position to id
X.id        = sparse([],[],[],M0,M1,X.maxcones) ;

% which ids are already assigned to cones
X.taken_ids = false(X.maxcones,1) ;

% four cardinal proximity relations, indexed by id
for d=1:4
    X.contact{d} = sparse(false(X.maxcones)) ;
end

end