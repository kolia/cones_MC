function X = initialize_X(M0,M1,N_colors,cell_consts,STA_W,coneConv,colorDot,D,maxcones,beta)

X.masks     = make_masks(D) ;
X.M0        = M0 ;
X.M1        = M1 ;
X.N_colors  = N_colors ;
X.N_cones   = 0  ;
X.ll        = 0 ;
X.n_moves   = 0 ;

X.beta      = beta ;
X.cell_consts = cell_consts ;
X.coneConv  = coneConv ;
X.colorDot  = colorDot ;

% upper bound on anticipated number of cones
X.maxcones  = maxcones ;

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
    X.contact{d} = sparse(false(X.maxcones)) ;
end

% transitive closure of contact forces X.contact
for d=1:4
    X.reach{d} = sparse(false(X.maxcones)) ;
end

% FROM old flip_color_LL.m !!!
X.state     = sparse([],[],[],1,size(STA_W,2),maxcones) ;
X.invWW     = [] ;
X.overlaps  = [] ;
%     X.WW        = [] ;
X.positions = [] ;
X.colors    = [] ;
X.sumLconst = length(cell_consts) * log(2*pi) + sum(log(cell_consts)) ;
%     X = flip_color_LL( temp_X , find(X.state) , prior_ll , cell_consts , ...
%                       STA_W , coneConv , colorDot , sizes , beta ) ;

end