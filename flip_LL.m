function X = flip_LL( X , flips , PROB , T )
% X = flip_LL( X , flips , PROB , T )
%
% pardon my appearance, i've been optimized for speed, not prettiness
%
% Apply flips to configuration X, and update log-likelihood of X.
% Bits in X.state are flipped, and the inverse X.invWW of WW is updated.
% Some other book-keeping variables are stored in X, for speed.
% The matrix inverse is not recalculated each time:
% block matrix inverse update formulas are used to update X.invWW 
% incrementally, for speed.

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

% WC = PROB.cone_params.weight_C ;  % covariance of prior weights
% Wm = PROB.cone_params.weight_m ;  % mean of prior weights

for i=1:size(flips,1)
    
    x = flips(i,1) ;
    y = flips(i,2) ;
    c = flips(i,3) ;
    
    [posX,posY,colors] = find(X.state) ;

    if ~c && ~X.state(x,y)
        error('deleting nonexistent cone...')
    else        % cone deletion
        k = x + M0*(y-1) ;
    end
    j = sum( posX + M0*(posY-1) <= k ) ;
    
    % block matrix inverse update
    if ~c       % update inverse by deleting jth row/column
        inds       = [1:j-1 j+1:X.N_cones] ;
        X.N_cones  = X.N_cones - 1 ;
        
        invWW      = X.invWW(inds,inds) - ...
                     X.invWW(inds,j)*X.invWW(j,inds)/X.invWW(j,j) ;
        X.invWW    = invWW ;
        
        X.STA_W_state = X.STA_W_state(:, inds ) ;
        
        X.state(x,y)= 0 ;
        
        X.diff = [X.diff ; x y 0] ;

    else        % update inverse by adding row/column
        j           = j + 1 ;
        X.N_cones   = X.N_cones + 1 ;
        inds        = [1:j-1 j+1:X.N_cones] ;
        
        Wkinds  = [posX'-x ; posY'-y] ;
        ssx     = 1+mod(x-1,PROB.SS) ;
        ssy     = 1+mod(y-1,PROB.SS) ;
        
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
        
        Wkkc = PROB.coneConv(PROB.R+ssx,PROB.R+ssy,ssx,ssy) * PROB.colorDot(c,c) ;
        
%         if max(Wkstate) > Wkkc
%             error('This is absurd.')
%         end
        
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
        
        STA_W_state = X.STA_W_state ;
        X.STA_W_state = zeros(PROB.N_GC,X.N_cones) ;
        if ~isempty(inds)
            X.STA_W_state(:,inds) = STA_W_state ;
        end
        X.STA_W_state(:,j) = PROB.make_STA_W( k, c, PROB.STA, PROB.cone_params.colors ) ;

        X.state(x,y)       = c ;        
        X.diff = [X.diff ; x y c] ;
    end
end

% recalculate data log-likelihood
ll = calculate_LL( X , PROB , T ) ;
X.T   = T  ;
X.ll  = ll ;

end