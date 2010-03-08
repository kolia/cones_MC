function X = flip_conv_LL( X , flips , prior_ll , cell_consts , STA_W , coneConv , sizes , beta )
% X = flip_conv_LL( X , flips , prior_ll , cell_consts , STA_W , coneConv , sizes , beta )
% Apply flips to configuration X, and update log-likelihood of X.
% Bits in X.state are flipped, and the inverse X.invWW of WW is updated.
% The matrix inverse is not recalculated each time:
% block matrix inverse update formulas are used to update X.invWW 
% incrementally, for speed.

% default inverse temperature is 1
if nargin<8
    beta = 1 ;
end

Nconv = ceil(size(coneConv,1)/2) ;

% initialize inverse if necessary
if ~isfield(X,'invWW')
    temp_X.state = false(1,size(STA_W,2)) ;
    temp_X.invWW = [] ;
    X = flip_conv_LL( temp_X , find(X.state) , prior_ll , cell_consts , STA_W , coneConv , sizes , beta ) ;
end

% initialize best 10 configurations encountered to zero
if ~isfield(X,'best')
    X.best = zeros(10,size(STA_W,2)+1) ;
    X.best(1:10) = -Inf ;
end

% apply all the bit flips to X
for i=1:length(flips)
    k = flips(i) ;
    n = size(X.invWW,1) ;
    j = cumsum(X.state) ;
    
    % block matrix inverse update
    if X.state(k)     % update inverse by deleting jth row/column
        j = j(k) ;
        inds = [1:j-1 j+1:n] ;
        X.invWW = X.invWW(inds,inds) - ...
                  X.invWW(inds,j   )*X.invWW(j,inds)/X.invWW(j,j) ;
    else                % update inverse by adding row/column
        j = j(k) + 1 ;
        inds = [1:j-1 j+1:n+1] ;
        A = X.invWW ;
        
        % reconstruct W(k,X.state) using coneConv
        Wkk     = coneConv(Nconv,Nconv) ;
        [sI,sJ] = ind2sub(sizes,find(X.state)) ;
        [kI,kJ] = ind2sub(sizes,k) ;
        Wkinds  = [sI-kI ; sJ-kJ] ;
        Wkstate = zeros(1,sum(X.state)) ;
        where   = max(abs(Wkinds),[],1)<Nconv ;
        if sum(where)
            here = sub2ind(size(coneConv),Nconv+Wkinds(1,where),Nconv+Wkinds(2,where)) ;
            Wkstate(where) = coneConv(here) ;
        end
        
        r = Wkstate * A ;
        q = 1/( Wkk - r*Wkstate' ) ;
        c = A * Wkstate' * q ;
        X.invWW = zeros(n+1) ;
        X.invWW(inds,inds) = A+c*r ;
        X.invWW(inds,j)    = -c    ;
        X.invWW(j,inds)    = -r*q  ;
        X.invWW(j,j)       = q     ;

    end
    X.state(flips(i)) = ~X.state(flips(i)) ;
end
X.state = logical(X.state) ;

% recalculate data log-likelihood
if sum(X.state)>0
    X.data_ll = ( - length(cell_consts) * log( det(2.*pi.*X.invWW) ) + ...
                sum( cell_consts' .* sum( (STA_W(:,X.state) * X.invWW) .* STA_W(:,X.state) ,2) ) ...
                )/2 ;
else
    X.data_ll = 0 ;
end

% update log-likelihood
X.ll = beta * (X.data_ll + prior_ll(X)) ;

% update best 10 configurations encountered so far
i = 1 ;
while 1    
    if i>10
        X.best = [ X.best(2:end,:) ; X.ll X.state ] ;
        break
    elseif sum(X.state ~= X.best(i,2:end))>0
        if X.ll > X.best(i)
            i = i + 1 ;
        else
            X.best = [X.best(1:i-1,:) ; X.ll X.state ; X.best(i:10,:)] ;
            X.best = X.best(2:end,:) ;
            break
        end
    else
        break
    end
end

end