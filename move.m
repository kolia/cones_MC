function sample = move( X , PROB , T )
% sample = move( X , PROB )
% Draw a trial flips starting from state X; calculated forward and backward
% proposal probabilities.
% Two types of move are produced: additions of cones, and moves or changes
% in color of existing cones.
% First, a spatial location is drawn. It is drawn from existing cone
% locations with probability PROB.q, and from the set of locations that are
% available for adding a cone with probability (1-PROB.q).
% For the trial location, if no cone of any color is at that position,
% then the trial is a placement of a single randomly colored cone at that
% location.
% If a cone of any color is at the location, then the corresponding
% trial is a change of color or deletion with probability 1/7, or a move of
% the cone to a neighboring location, each with probability 1/7. If the
% cone lies on the border of the domain, 7 is replaced by 3 + the number of
% adjacent positions that are within the border.

if nargin<4 ,  n = 1 ; end

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

% current cones in X
[cx,cy]     = find(X.state) ;

% probability of changing existing cone
if X.N_cones >= X.maxcones  % not adding any more cones
    q = 1 ;
else
    q = PROB.q ;
end


% propose moves of existing cones
if X.N_cones>0  &&  rand()<q
    
    cone = randi(X.N_cones,1) ;
    
    i       = cx(cone) ;
    j       = cy(cone) ;
    color   = X.state(i,j) ;

    % number of legal moves for this cone, being careful about borders
    nforward = 4 - PROB.outofbounds(i,j) + PROB.N_colors ;

    % cone deletion
    if rand()<1/nforward                    
        sample = flip_LL( X , [i j 0] , PROB , T ) ;

        % number of available new cone locations after removing cone
        [~,indices] = not_excluded( X, i, j ) ;
        excluded = X.excluded ;
        excluded(indices) = false ;
        N_not_excluded = M0*M1 - nnz(excluded) ;
        
        forward_probability = 1/X.N_cones * q / nforward ;

        % probability of proposing reverse move: re-adding cone
        reverse_probability = (1-PROB.q)/(N_not_excluded*PROB.N_colors) ;
        
        % proposal probability bias
        sample.proposal_bias = forward_probability/reverse_probability ;

    else
        % change cone color
        if rand()<(PROB.N_colors-1)/(nforward-1)
            new_color = setdiff( 1:X.N_colors , color ) ;
            new_color = new_color( randi(numel(new_color)) ) ;

            ne = not_excluded(X,i,j) ;
            if ne && ~isempty(PROB.sparse_struct{i,j,new_color})
                sample = flip_LL( X , [i j new_color] , PROB , T ) ;
            else
                sample    = X ;
                sample.ll = -Inf ;
            end
            
        % shift cone    
        else
            
            % choose shift direction that is within bounds
            while 1
                d = randi(4,1) ;
                ni = i + X.masks{1,1}.shift{d}(1) ;
                nj = j + X.masks{1,1}.shift{d}(2) ;
                if ni > 0   &&  ni <= M0  &&  nj>0  &&  nj <= M1
                    break ;
                end
            end
            sample = propagate_action(X,i+(j-1)*X.M0,d,PROB,T) ;
        end
        
        % proposal probability bias
        sample.proposal_bias = 1 ;
    end
    
else
    
    % sample unoccupied locations, propose cone additions
    possible = find(~X.excluded & PROB.has_evidence) ;
    [i,j] = ind2sub([M0 M1],possible(randi(numel(possible)))) ;

    % forward probability of choosing this location & color
    forward_probability = (1-q)/(numel(possible) * PROB.N_colors) ;

    % choose new color
    new_color = randi(PROB.N_colors) ;

    % propose addition of new cone of each color
    ne = not_excluded(X,i,j) ;
    if ne && ~isempty(PROB.sparse_struct{i,j,new_color})
        sample = flip_LL( X , [i j new_color] , PROB , T ) ;
    else
        sample    = X ;
        sample.ll = -Inf ;                
    end
    
    % probability of proposing reverse move: re-deleting cone
    if X.N_cones >= X.maxcones  % not adding any more cones
        q = 1 ;
    else
        q = PROB.q ;
    end
    nforward = 4 - PROB.outofbounds(i,j) + PROB.N_colors ;
    reverse_probability = 1/X.N_cones * q / nforward ;
    
    % proposal probability bias
    sample.proposal_bias = forward_probability/reverse_probability ;
end

end