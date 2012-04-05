function samples = move( X , PROB , T , n )
% samples = move( X , PROB , n )
% Draw n trial flips starting from state X.
% Two types of move are produced: additions of cones, and moves or changes
% in color of existing cones.
% First, n spatial locations are drawn. They are drawn from existing cone
% locations with probability PROB.q, and from the uniform distribution with
% probability (1-PROB.q).
% For each trial location i, if no cone of any color is at that position,
% then the corresponding trial is a placement of a single randomly colored
% cone at that location.
% If a cone of any color is at location i, then the corresponding
% trial is a change of color or deletion with probability 1/7, or a move of
% the cone to a neighboring location, each with probability 1/7. If the
% cone lies on the border of the domain, 7 is replaced by 3 + the number of
% adjacent positions that are within the border.
% For speed, backward probabilities are not calculated: this sampler should
% only be used with the symmetric rule, not with metropolis-hastings.

% check_X(X)

if nargin<4 ,  n = 1 ; end

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

% current cones in X
[cx,cy]     = find(X.state) ;

% initialize samples
samples  = cell(X.N_cones*7+n,1) ;
ns       = 0 ;

% propose moves of existing cones
if X.N_cones > 0
    % draw n_moved existing cones
    n_moved     = binornd(n,PROB.q) ;
    cones       = randi( X.N_cones , 1 , n_moved ) ;
    
%     cones       = 1:X.n_cones ;
    
    for s=1:n_moved
%     for s=1:X.n_cones
        i       = cx(cones(s)) ;
        j       = cy(cones(s)) ;
        color   = X.state(i,j) ;
        
        % probability of choosing this location
        p = 1/X.N_cones * PROB.q ;
        
        % number of legal moves for this cone, being careful with borders
        nforward    = 4 - PROB.outofbounds(i,j) + PROB.N_colors ;

        % for each adjacent location, add trial move to that location
        for d=1:4                 % move N , E , S , W
            ni = i + X.masks{1,1}.shift{d}(1) ;
            nj = j + X.masks{1,1}.shift{d}(2) ;
            if ni > 0   &&  ni <= M0  &&  nj>0  &&  nj <= M1
                ns = ns+1 ;
%                 samples{ns} = move_cone( X , i , j , d , PROB , T ) ;                
                samples{ns} = propagate_action(X,i+(j-1)*X.M0,d,PROB,T) ;
                samples{ns}.forward_prob   = p/nforward ;
            end
        end
        
        % cone deletion
        ns = ns+1 ;
        samples{ns} = flip_LL( X , [i j 0] , PROB , T ) ;
        samples{ns}.forward_prob    = p/nforward ;
%         check_X(samples{ns})

        % change of color, without moving
        deleted = ns ;
        ne = not_excluded(X,i,j) ;
        for cc=setdiff( 1:X.N_colors , color )
            ns = ns+1 ;
            if ne && ~isempty(PROB.sparse_struct{i,j,cc})
                samples{ns} = flip_LL( samples{deleted} , [i j cc] , PROB , T ) ;
%                 check_X(samples{ns})
            else
                samples{ns}    = samples{deleted} ;
                samples{ns}.ll = -Inf ;
            end
            samples{ns}.forward_prob    = p/nforward ;
        end
        
    end
else
    n_moved = 0 ;
end

% sample unoccupied locations, propose cone additions
while ns <= n - n_moved
    ex = find(~X.excluded & PROB.has_evidence) ;
    [i,j] = ind2sub([M0 M1],ex(randi(numel(ex)))) ;
%     i = ij(1) ;
%     j = ij(2) ;
%     i       = randi( M0 , 1 ) ;
%     j       = randi( M1 , 1 ) ;
    
    if ~X.state(i+M0*(j-1))
        % probability of choosing this location & color
        p = (1-PROB.q)/((M0*M1 - X.N_cones)*PROB.N_colors) ;

        if X.N_cones >= X.maxcones
        % if maxcones has been reached, delete a random cone first
            cone = randi( X.N_cones , 1 , 1 ) ;
            X    = flip_LL( X, [cx(cone) cy(cone) 0], PROB, T ) ;
            p    = p/X.N_cones ;
        end
        
        % propose addition of new cone of each color
        ne = not_excluded(X,i,j) ;
        for c=1:PROB.N_colors
            ns = ns+1 ;
            if ne && ~isempty(PROB.sparse_struct{i,j,c})
                samples{ns} = flip_LL( X , [i j c] , PROB , T ) ;
%                 check_X(samples{ns})
            else
                samples{ns}    = X ;
                samples{ns}.ll = -Inf ;                
            end
            samples{ns}.forward_prob    = p ;
        end
    end
end

samples = samples(1:ns) ;
% fprintf('\nDONE sampling trial moves\n')

end