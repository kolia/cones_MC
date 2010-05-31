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
    x = flips(i,1) ;
    y = flips(i,2) ;
    c = flips(i,3) ;
    
    j = sum( find(X.state) <= x + X.M0*(y-1) + X.M0*X.M1*(c-1) ) ;
    
    % block matrix inverse update
    if X.state(k)     % update inverse by deleting jth row/column
        inds = [1:j-1 j+1:X.N_cones] ;
        X.overlaps = X.overlaps(inds,inds) ;
        X.invWW = X.invWW(inds,inds) - X.invWW(inds,j)*X.invWW(j,inds)/X.invWW(j,j) ;
        %         X.WW    = X.WW(inds,inds) ;
        X.positions = X.positions(:,inds) ;
        X.colors    = X.colors(inds) ;
        X.state(k)  = false ;
        
    else               % update inverse by adding row/column
        j = j + 1 ;
        inds = [1:j-1 j+1:X.N_cones+1] ;
        A = X.invWW ;
        
        % reconstruct W(k,X.state) using coneConv
        Wkk     = coneConv(Nconv,Nconv) ;
        
        positions   = X.positions ;
        X.positions = zeros(2,X.N_cones + 1) ;
        if ~isempty(positions)
            X.positions(:,inds) = positions ;
            Wkinds  = [positions(1,:)-x ; positions(2,:)-y] ;
        else
            Wkinds  = [] ;
        end
        X.positions(:,j)    = [x;y] ;
        
        colors      = X.colors ;
        X.colors    = zeros(1,X.N_cones + 1) ;
        X.colors(inds) = colors ;
        X.colors(j)    = c  ;
        
        overlap = zeros(1,X.N_cones) ;
        Wkstate = overlap ;
        where   = max(abs(Wkinds),[],1)<Nconv ;
        if sum(where)
            Nc   = size(coneConv,1) ;
            here = (Nconv+Wkinds(2,where) - 1)*Nc + Nconv+Wkinds(1,where) ;
            overlap(where) = coneConv(here) ;
            Wkstate(where) = overlap(where) .* colorDot(c,colors(where)) ;
        end
        
        Wkkc = Wkk * colorDot(c,c) ;
        
        %         WW                  = X.WW ;
        %         X.WW                = zeros(n+1) ;
        %         X.WW(inds,inds)     = WW ;
        %         X.WW(inds,j)        = Wkstate ;
        %         X.WW(j,inds)        = Wkstate ;
        %         X.WW(j,j)           = Wkkc ;
        
        O                     = X.overlaps ;
        X.overlaps            = zeros(X.N_cones+1) ;
        X.overlaps(inds,inds) = O ;
        X.overlaps(inds,j)    = overlap ;
        X.overlaps(j,inds)    = overlap ;
        X.overlaps(j,j)       = Wkk ;
        
        r = Wkstate * A ;
        q = 1/( Wkkc - r*Wkstate' ) ;
        c = A * Wkstate' * q ;
        X.invWW = zeros(X.N_cones+1) ;
        X.invWW(inds,inds) = A+c*r ;
        X.invWW(inds,j)    = -c    ;
        X.invWW(j,inds)    = -r*q  ;
        X.invWW(j,j)       = q     ;
        
        X.state(k)  = true ;
        
    end
end
X.state = logical(X.state) ;

% recalculate data log-likelihood
Ncones = sum(X.state) ;
if Ncones>0
    
    invWW = X.invWW ;
    invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;
    invWW = sparse(invWW) ;
    ldet  = 2 * sum(log(diag(chol(invWW))));
    
    STA_W_state = STA_W(:,X.state) ;

    X.ll  = X.beta * full(- Ncones * X.sumLconst + ldet + ...
        sum( cell_consts .* sum( (STA_W_state * invWW) .* STA_W_state ,2) )/2) ;

else
    X.ll = 0 ;
end

end