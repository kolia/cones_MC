function X = flip_LL( X , flips , cell_consts , STA_W )
% X = flip_LL( X , flips , cell_consts , STA_W )
%
% pardon my appearance, i've been optimized for speed, not prettiness
%
% Apply flips to configuration X, and update log-likelihood of X.
% Bits in X.state are flipped, and the inverse X.invWW of WW is updated.
% Some other book-keeping variables are stored in X, for speed.
% The matrix inverse is not recalculated each time:
% block matrix inverse update formulas are used to update X.invWW 
% incrementally, for speed.

% center of cone RF convolution
Nconv   = ceil(size(X.coneConv,1)/2) ;

for i=1:size(flips,1)

%     OX = X ;
    
    x = flips(i,1) ;
    y = flips(i,2) ;
    c = flips(i,3) ;
    
    [posX,posY,colors] = find(X.state) ;

    if ~c && ~X.state(x,y)
        error('deleting nonexistent cone...')
    elseif c    % cone addition
        k = x + X.M0*(y-1) + X.M0*X.M1*(c-1) ;
    else        % cone deletion
        k = x + X.M0*(y-1) + X.M0*X.M1*(X.state(x,y)-1) ;
    end
    j = sum( posX + X.M0*(posY-1) + X.M0*X.M1*(colors-1) <= k ) ;
    
    % block matrix inverse update
    if ~c       % update inverse by deleting jth row/column
        inds = [1:j-1 j+1:X.N_cones] ;
        X.overlaps = X.overlaps(inds,inds) ;
        
        invWW      = X.invWW(inds,inds) - X.invWW(inds,j)*X.invWW(j,inds)/X.invWW(j,j) ;
        if c
            X.invWW(inds,inds) = invWW ;
        else
            X.invWW = invWW ;
        end
        %         X.WW    = X.WW(inds,inds) ;
        X.state(x,y)= 0 ;
        X.N_cones   = X.N_cones - 1 ;
        
        X.diff.deleted = [X.diff.deleted ; x y] ;
        
    else                     % update inverse by adding row/column
        j = j + 1 ;
        X.N_cones     = X.N_cones + 1 ;
        inds = [1:j-1 j+1:X.N_cones] ;
        A = X.invWW ;
                
        % reconstruct W(k,X.state) using coneConv
        Wkk     = X.coneConv(Nconv,Nconv) ;
        
        Wkinds  = [posX-x ; posY-y] ;
        
        overlap = zeros(1,X.N_cones-1) ;
        Wkstate = overlap ;
        where   = max(abs(Wkinds),[],1)<Nconv ;
        if sum(where)
            Nc   = size(X.coneConv,1) ;
            here = (Nconv+Wkinds(2,where) - 1)*Nc + Nconv+Wkinds(1,where) ;
            overlap(where) = X.coneConv(here) ;
            Wkstate(where) = overlap(where) .* X.colorDot(c,colors(where)) ;
        end
        
        Wkkc = Wkk * X.colorDot(c,c) ;
        
        %         WW                  = X.WW ;
        %         X.WW                = zeros(n+1) ;
        %         X.WW(inds,inds)     = WW ;
        %         X.WW(inds,j)        = Wkstate ;
        %         X.WW(j,inds)        = Wkstate ;
        %         X.WW(j,j)           = Wkkc ;
        
        O                     = X.overlaps ;
        X.overlaps            = zeros(X.N_cones) ;
        X.overlaps(inds,inds) = O ;
        X.overlaps(inds,j)    = overlap ;
        X.overlaps(j,inds)    = overlap ;
        X.overlaps(j,j)       = Wkk ;
        
        r = Wkstate * A ;
        
        X.invWW = zeros(X.N_cones) ;
        if ~isempty(r)
            q = 1/( Wkkc - r*Wkstate' ) ;
            cc= A * Wkstate' * q ;
            X.invWW(inds,inds) = A+cc*r  ;
            X.invWW(inds,j)    = -cc     ;
            X.invWW(j,inds)    = -r*q    ;
        else
            q = 1/Wkkc ;
        end
        X.invWW(j,j)       = q ;
        X.state(x,y)       = c ;
        
        X.diff.added = [X.diff.added ; x y] ;

    end
end

% recalculate data log-likelihood
if X.N_cones>0
    
    invWW = X.invWW ;
    invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;
    invWW = sparse(invWW) ;
    ldet  = 2 * sum(log(diag(chol(invWW))));
    
    [x,y,c] = find(X.state) ;
    
    STA_W_state = STA_W( x+X.M0*(y-1)+X.M0*X.M1*(c-1) , : )' ;
    
%     try
    ll  = X.beta * full(- X.N_cones * X.sumLconst + ldet + ...
        sum( cell_consts .* sum( (STA_W_state * invWW) .* STA_W_state ,2) )/2) ;
%     catch
%        'blah' 
%     end
else
    ll = 0 ;
end

X.dll = ll - X.ll ;
X.ll  = ll ;

end