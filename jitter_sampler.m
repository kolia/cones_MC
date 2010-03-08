function samples = jitter_sampler( X , p , n , M )
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

sampled = randiscrete( p , n) ;
samples = cell(n*8,1) ;
ns      = 0 ;

for s=1:n
    k  = sampled(s) ;
    if X.state(k) % 
        [i,j]= ind2sub([M M],k) ;
        sx   = max(1,i-1):min(M,i+1) ;
        sy   = max(1,j-1):min(M,j+1) ;
        sxy  = [reshape( repmat( sx  , length(sy) , 1 ) , [] , 1 ) ...
            reshape( repmat( sy' , 1 , length(sx) ) , [] , 1 ) ] ;
        inds = sub2ind( [M M] , sxy(:,1) , sxy(:,2) ) ;
                
        ni = length(inds) ;
        for i=1:ni
            ns = ns+1 ;
            if inds(i) ~= k
                samples{ns}.flips = [k inds(i)] ;
            else
                samples{ns}.flips = k ;
            end
            samples{ns}.forward_prob   = p(k)/ni ;

%             nb  = (mod( sxy(i,1) , M ) < 2)  +  (mod( sxy(i,2) , M ) < 2) ;
%             key = [9 6 4] ;
%             samples{ns}.backward_prob  = p(k)/key(nb+1) ;
        end
    else
        ns = ns+1 ;
        samples{ns}.flips = k ;
        samples{ns}.forward_prob  = p(k) ;
        
%         nb  = (mod( sxy(5,1) , M ) < 2)  +  (mod( sxy(5,2) , M ) < 2) ;
%         key = [9 6 4] ;
%         samples{ns}.backward_prob = p(k)/key(nb+1) ;
    end
end

samples = samples(1:ns) ;