function X = initialize_X(M0,M1,N_colors,SS,repulsion,naive_LL,beta,delta,maxcones)

if nargin<9
    maxcones = floor( 0.3 * M0 * M1 ) ;
end

X.masks     = make_masks(repulsion*SS) ;
X.D         = max(repulsion(1:2,1)*SS) ;

X.repulsion = repulsion ;
X.M0        = M0 * SS ;
X.M1        = M1 * SS ;
X.SS        = SS ;
X.N_colors  = N_colors ;
X.N_cones   = 0  ;
X.ll        = 0 ;
X.dll       = 0 ;
X.n_moves   = 0 ;
X.diff      = [] ;
X.version   = 0 ;
X.naive_LL  = naive_LL ;
X.excluded  = false(X.M0,X.M1) ;

X.beta      = beta  ;
X.delta     = delta ;

% upper bound on anticipated number of cones
X.maxcones  = maxcones ;

% sparse int matrix, representing cone positions and colors
X.state     = sparse([],[],[],X.M0,X.M1,X.maxcones) ;

% X.STA_W_state = [] ;
X.sparse_STA_W_state = sparse([]) ;

% contact forces at four cardinal adjacent positions, indexed by id
for d=1:2
    X.contact{d} = sparse([],[],[],X.M0*X.M1,X.M0*X.M1,floor(X.maxcones/2)) ;
    X.contact{d} = logical(X.contact{d}) ;
end

nn = ceil(M0*M1/10) ;
X.cputime = zeros(nn,1) ;
X.LL_history= zeros(nn,1) ;
X.N_cones_history = zeros(nn,1) ;
X.T = [1 1] ;

% FROM old flip_color_LL.m !!!
% X.invWW     = [] ;
% X.overlaps  = [] ;
X.WW        = [] ;
% X.positions = [] ;
% X.colors    = [] ;

end