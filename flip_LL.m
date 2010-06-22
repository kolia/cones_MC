function X = flip_LL( X , flips , PROB )
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

for i=1:size(flips,1)
    
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
        inds       = [1:j-1 j+1:X.N_cones] ;
        X.N_cones  = X.N_cones - 1 ;
        
        invWW      = X.invWW(inds,inds) - X.invWW(inds,j)*X.invWW(j,inds)/X.invWW(j,j) ;
        X.invWW    = invWW ;
        
        X.state(x,y)= 0 ;
        
        X.diff.deleted = [X.diff.deleted ; x y] ;
        
    else        % update inverse by adding row/column
        j           = j + 1 ;
        X.N_cones   = X.N_cones + 1 ;
        inds        = [1:j-1 j+1:X.N_cones] ;
        
        Wkinds  = [posX'-x ; posY'-y] ;
        ssx     = 1+mod(1-x,X.SS) ;
        ssy     = 1+mod(1-y,X.SS) ;
        
%         fprintf('\t%d',min(max(abs(Wkinds),[],1)))
        
        Wkstate = zeros(1,X.N_cones-1) ;
        where   = find( max(abs(Wkinds),[],1) <= PROB.R ) ;
        if ~isempty(where)
            xx = Wkinds(1,where)+PROB.R+ssx ;
            yy = Wkinds(2,where)+PROB.R+ssy ;
            for k=1:length(where)
                Wkstate(where(k)) = PROB.coneConv(xx(k),yy(k),ssx,ssy) ...
                    .* PROB.colorDot(c,colors(where(k))) ;
            end
        end
        
        Wkkc    = PROB.coneConv(PROB.R+1,PROB.R+1,ssx,ssy) * PROB.colorDot(c,c) ;        
        
        if max(Wkstate) > Wkkc
            error('This is absurd.')
        end
        
        invWW   = X.invWW ;
        r       = Wkstate * invWW ;
        
        X.invWW = zeros(X.N_cones) ;
        if ~isempty(r)
            q                  = 1/( Wkkc - r*Wkstate' ) ;
            cc                 = invWW * Wkstate' * q ;
            X.invWW(inds,inds) = invWW+cc*r  ;
            X.invWW(inds,j)    = -cc         ;
            X.invWW(j,inds)    = -r*q        ;
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
    
    STA_W_state = PROB.STA_W( x+X.M0*(y-1)+X.M0*X.M1*(c-1) , : )' ;
    
    ll  = X.beta * full(- X.N_cones * PROB.sumLconst + ldet * PROB.N_GC + ...
        sum( PROB.cell_consts .* sum( (STA_W_state * invWW) .* STA_W_state ,2) )/2) ;
else
    ll = 0 ;
end

X.dll = ll - X.ll ;
X.ll  = ll ;

end