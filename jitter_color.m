function samples = jitter_color( X , n , q , flip_LL , N_colors )
% samples = jitter_color( X , n , q , flip_LL , N_colors )
% Draw n trial flips starting from state X.
% M0 -by- M1 is the spatial extent of cone locations, N-colors = 3. 
% Two types of move are produced: additions of cones, and moves or changes
% in color of existing cones.
% First, n spatial locations are drawn. They are drawn from existing cone
% locations with probability q, and from the uniform distribution with
% probability (1-q).
% For each trial location i, if no cone of any color is at that position,
% then the corresponding trial is a placement of a single randomly colored
% cone at that location, or adding no cones at all (nothing changes).
% If a cone of any color is at location i, then the corresponding 
% trial is a change of color or deletion with probability 1/9, or a move of
% the cone to a neighboring location, each with probability 1/9. If the
% cone lies on the border of the domain, 9 is replaced by 1 + the number of
% available adjacent bits.
% For speed, backward probabilities are not calculated: this sampler should
% only be used with the symmetric rule, not with metropolis-hastings.

% defaults
if nargin<7 ,   N_colors = 3 ;                end

[M0,M1] = size(X.state) ;

% current cones in X
[cx,cy]     = find(X.state>0) ;
n_cones     = numel( cx ) ;

% initialize samples
samples  = cell(n*8,1) ;
ns       = 0 ;

% number of draws from existing cone locations
n_moved = binornd(n,q) ;

% draw cone locations
if n_cones > 0
    cones       = randi( n_cones , 1 , n_moved ) ;

    for s=1:n_moved
        i       = cx(cones(s)) ;
        j       = cy(cones(s)) ;
        color   = X.state(i,j) ;
        
        % probability of choosing this location
        p = 1/n_cones * q ;

        % adjacent locations, being careful not to go beyond borders
        sx   = max(1,i-1):min(M0,i+1) ;
        sy   = max(1,j-1):min(M1,j+1) ;

        % for each adjacent location, add trial move to that location
        ni = numel(sx) * numel(sy) ;
        Xtemp = flip_LL(X,i,j,color) ;
        for x=sx
            for y=sy
                % move to adjacent location
                if x~=i || y~=j
                    ns = ns+1 ;
                    samples{ns} = flip_LL(Xtemp,x,y,color) ;
                    samples{ns}.forward_prob   = p/ni ;
                    
                % change of color or deletion, without moving
                else
                    for cc=1:N_colors
                        ns = ns+1 ;
                        samples{ns} = flip_LL(X,x,y,cc) ;
                        samples{ns}.forward_prob = p/(ni*N_colors) ;
                    end
                end
            end
        end
    end
else
    n_moved = 0 ;
end


sampled_x = randi( M0 , n-n_moved) ;
sampled_y = randi( M1 , n-n_moved) ;

% for each sampled location, generate corresponding moves
for s=1:n-n_moved
    i       = sampled_x(s) ;
    j       = sampled_y(s) ;
    
    % probability of choosing this location
    p = (1-q)/(M0*M1) ;
    
    % store away change of color
    for cc=1:N_colors
        ns = ns+1 ;
        samples{ns} = flip_LL(X,i,j,cc) ;
        samples{ns}.forward_prob  = p/N_colors ;
    end
end

% for s=1:ns
%     k = samples{s}.flips(end) ;
%     position = mod( k-1 , NBW ) + 1 ;
%     [i,j]    = ind2sub( [M0 M1] , position ) ;
%     color    = floor((k-1) / NBW) + 1 ;
%     fprintf('\nPosition %3d %3d  Color %d',i,j,color)
% end
% fprintf('\n')

samples = samples(1:ns) ;

end