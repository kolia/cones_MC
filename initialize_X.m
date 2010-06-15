function X = initialize_X(M0,M1,N_colors,SS,maxcones,D,beta)

X.masks     = make_masks(D) ;
X.M0        = M0 * SS ;
X.M1        = M1 * SS ;
X.SS        = SS ;
X.N_colors  = N_colors ;
X.N_cones   = 0  ;
X.ll        = 0 ;
X.dll       = 0 ;
X.n_moves   = 0 ;
X.diff.added= [] ;
X.diff.deleted = [] ;

X.beta      = beta ;

% upper bound on anticipated number of cones
X.maxcones  = maxcones ;

% sparse int matrix, representing cone positions and colors
X.state     = sparse([],[],[],X.M0,X.M1,X.maxcones) ;

% sparse int matrix, with number of out-of-border adjacencies
X.outofbounds = sparse([],[],[],X.M0,X.M1,2*(X.M0+X.M1)) ;
X.outofbounds(:,[1 X.M1]) = 1 ;
X.outofbounds([1 X.M0],:) = X.outofbounds([1 X.M0],:) + 1 ;

% contact forces at four cardinal adjacent positions, indexed by id
for d=1:4
    X.contact{d} = sparse([],[],[],X.M0*X.M1,X.M0*X.M1,X.maxcones) ;
    X.contact{d} = logical(X.contact{d}) ;
end

% FROM old flip_color_LL.m !!!
X.invWW     = [] ;
X.overlaps  = [] ;
%     X.WW        = [] ;
X.positions = [] ;
X.colors    = [] ;

end