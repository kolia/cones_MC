function [X,swaps] = swap_randrect( X , otherX , n , XLL , OLL )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations:  swaps of random rectangular regions of space are stored
% in swap.flips. 

[M0,M1,N_colors] = size(X.state) ;
swaps = cell(n,1) ;

ns = 0 ;
while ns<n
    
    xrange = randrange(M0,40) ;
    yrange = randrange(M1,40) ;
    
    Xcones =      X.state(xrange(1):xrange(2),yrange(1):yrange(2)) ;
    Ocones = otherX.state(xrange(1):xrange(2),yrange(1):yrange(2)) ;
    
    [xsupport,ysupport] = find( xor( Xcones>0 , Ocones>0 ) ) ;

    if ~isempty(xsupport)
        ns = ns + 1 ;

        YX = X ;
        YO = otherX ;
        for i=1:length(xsupport)
            x = xrange(1) + xsupport(i) - 1 ;
            y = yrange(1) + ysupport(i) - 1 ;
            YX = XLL(YX,x,y,otherX.state(x,y)) ;
            YO = OLL(YO,x,y,     X.state(x,y)) ;
        end
        swaps{ns}.X     = YX ;
        swaps{ns}.with  = YO ;
        swaps{ns}.ll    = YX.ll + YO.ll ;
        swaps{ns}.state = YX.state ;
        swaps{ns}.forward_prob = 1 ;  % all swaps have the same prob anyway
    end
end

X.X     = X ;
X.with  = otherX ;
X.ll    = X.X.ll + X.with.ll ;
X.state = X.X.state ;

end


function range = randrange(M,m)

n = max(2,ceil(M/m)) ;
r = randi(M,1,n) ;
r = sort(r) ;
i = randi(n-1) ;
range = [r(i) r(i+1)] ;

end