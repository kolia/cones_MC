function samples = jitter_color( X , n , q , M0 , M1 , cumprob , N_colors )
% samples = jitter_color( X , n , q , M0 , M1 , cumprob , N_colors )
% Draw n trial flips starting from state X.
% M0 -by- M1 is the spatial extent of cone locations, N-colors = 3. 
% Two types of move are produced: additions of cones, and moves or changes
% in color of existing cones.
% First, n spatial locations are drawn. They are drawn from existing cone
% locations with probability q, and from probability cumprob with
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
if nargin<6 ,   cumprob  = (1:M0*M1) * 2 / (M0*M1*(M0*M1+1)) ;   end

% current number of cones in X
n_cones = size(X.positions,2) ;

% number of spatial locations, should equal  M0 * M1
NBW = numel(X.state)/N_colors ;

% number of draws from existing cone locations
n_moved = binornd(n,q) ;

% draw cone locations
sampled = randiscrete( cumprob , n-n_moved) ;
if n_cones > 0
    cones   = randi( n_cones , 1 , n_moved ) ;
    sampled  = [X.positions(1,cones)+M0*(X.positions(2,cones)-1) sampled] ;
else
    sampled = randiscrete( cumprob , n) ;
end

% initialize samples
samples  = cell(n*8,1) ;
ns       = 0 ;

% for each sampled location, generate corresponding moves
for s=1:length(sampled)
    
    % current position index considered
    k  = sampled(s) ;
    
    % compute corresponding i-j coordinates and color
    position = mod( k-1 , NBW      ) + 1 ;
    color    = mod( k-1 , N_colors ) + 1 ;
    i        = 1 + mod( position-1 , M0 ) ;
    j        = 1 + floor( (position-1) / M0) ;
    
    % does this position already have a cone?
    occupied= s<=n_moved  &&  ismember([i j],X.positions','rows') ;
    
    % probability of choosing this spatial location from cumprob
    if position == 1
        p = cumprob(1) ;
    else
        p = cumprob(position) - cumprob(position-1) ;
    end
    
    
    % if location already had a cone
    if occupied

        % probability of choosing this location
        p = 1/n_cones * q  +  (1-q) * p ;

        % adjacent locations, being careful not to go beyond borders
        sx   = max(1,i-1):min(M0,i+1) ;
        sy   = max(1,j-1):min(M1,j+1) ;
        sxy  = [reshape( repmat( sx  , length(sy) , 1 ) , [] , 1 ) ...
            reshape( repmat( sy' , 1 , length(sx) ) , [] , 1 ) ] ;
        inds = sub2ind( [M0 M1] , sxy(:,1) , sxy(:,2) ) ;

        % for each adjacent location, add trial move to that location
        ni = length(inds) ;
        for i=1:ni
            % store away move to adjacent location
            if inds(i) ~= position
                ns = ns+1 ;
                samples{ns}.flips = [position inds(i)] + (color-1)*M0*M1 ;
                samples{ns}.forward_prob   = p/ni ;

            % store away change of color or deletion, without moving
            else
                for j=1:N_colors
                    ns = ns+1 ;
                    samples{ns}.flips = position + (j-1)*M0*M1 ;
                    samples{ns}.forward_prob = p/ni/N_colors ;
                end
            end
        end
    
    % if location unoccupied by cone
    else
        % probability of choosing this location
        p = (1-q)*p ;

        % store away change of color
        for j=1:N_colors
            ns = ns+1 ;
            samples{ns}.flips = position + (j-1)*M0*M1 ;
            samples{ns}.forward_prob  = p/N_colors ;
        end
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