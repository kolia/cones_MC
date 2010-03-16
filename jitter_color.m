function samples = jitter_color( X , n , M0 , M1 , cumprob , N_colors )
% Draw n trial flips starting from state X.
% Two types of flips are produced: single bit flips, and shifts (moves).
% First, n samples from discrete probability p are drawn.
% For each sample i, if X.state(i) is 0, then the corresponding trial flip
% is a single flip of bit i. If X.state(i) is 1, then the corresponding
% trial flip is a single flip of bit i with probability 1/9, or a move of
% this 1 to a neighboring location, each with probability 1/9. If bit i
% lies on the border of the domain, 9 is replaced by the number of
% available adjacent bits.
% For speed, backward probabilities are not calculated: this sampler should
% only be used with the symmetric rule, not with metropolis-hastings.

if nargin<6 ,   N_colors = 3 ;                end
if nargin<5 ,   p = ones(M0*M1,1)/(M0*M1) ;   end

sampled  = randiscrete( cumprob , n) ;
samples  = cell(n*8,1) ;
ns       = 0 ;

% occupied = sum( reshape( X.state , [] , N_colors ) , 2) ;

NBW = numel(X.state)/N_colors ;

for s=1:n
    k  = sampled(s) ;
    position = mod( k-1 , NBW      ) + 1 ;
    color    = mod( k-1 , N_colors ) + 1 ;
    if position == 1
        p = cumprob(1) ;
    else
        p = cumprob(position) - cumprob(position-1) ;
    end

    if ismember(X.positions,position)
        [i,j]= ind2sub( [M0 M1] , position ) ;
        sx   = max(1,i-1):min(M0,i+1) ;
        sy   = max(1,j-1):min(M1,j+1) ;
        sxy  = [reshape( repmat( sx  , length(sy) , 1 ) , [] , 1 ) ...
            reshape( repmat( sy' , 1 , length(sx) ) , [] , 1 ) ] ;
        inds = sub2ind( [M0 M1] , sxy(:,1) , sxy(:,2) ) ;

        ni = length(inds) ;
        for i=1:ni
            if inds(i) ~= position
                ns = ns+1 ;
                samples{ns}.flips = [position inds(i)] + (color-1)*M0*M1 ;
                samples{ns}.forward_prob   = p/ni ;
            else
                for j=1:N_colors
                    ns = ns+1 ;
                    samples{ns}.flips = position + (j-1)*M0*M1 ;
                    samples{ns}.forward_prob = p/ni/N_colors ;
                end
            end
        end
    else
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