function X = initialize_X(M0,M1,N_colors,SS,repulsion,beta,delta,maxcones)

% a decent rule-of-thumb upper bound on the possible number of cones
if nargin<9,  maxcones = floor( 0.3 * M0 * M1 ) ; end

% masks is a data structure for calculating cone exclusion and contacts
X.masks     = make_masks(repulsion*SS) ;

% D is the exclusion distance used throughout
X.D         = max(repulsion(1:2,1)*SS) ;

X.repulsion = repulsion ;

% size of the ROI
X.M0        = M0 * SS ;
X.M1        = M1 * SS ;

% supersampling factor; usually 4
X.SS        = SS ;

% 3 colors
X.N_colors  = N_colors ;

% initial configuration
X.N_cones   = 0  ;
X.ll        = 0 ;
X.dll       = 0 ;
X.n_moves   = 0 ;
X.diff      = [] ;
X.version   = 0 ;
X.excluded  = false(X.M0,X.M1) ;

% temperature progression, used by CAST
X.beta      = beta  ;
X.delta     = delta ;

% upper bound on anticipated number of cones
X.maxcones  = maxcones ;

% cone positions and colors, the most important data structure in X
X.state     = sparse([],[],[],X.M0,X.M1,X.maxcones) ;

% used by flip_LL to calculate all log-likelihoods
X.WW        = [] ;
X.sparse_STA_W_state = sparse([]) ;

% horizontal and vertical cone contact relations, used during shift moves
for d=1:2
    X.contact{d} = sparse([],[],[],X.M0*X.M1,X.M0*X.M1,floor(X.maxcones/2)) ;
    X.contact{d} = logical(X.contact{d}) ;
end

% record cputime, LL_history, and N_cones_history after each iteration
nn = ceil(M0*M1/10) ;
X.cputime = zeros(nn,1) ;
X.LL_history= zeros(nn,1) ;
X.N_cones_history = zeros(nn,1) ;

% default temperature
X.T = [1 1] ;


end