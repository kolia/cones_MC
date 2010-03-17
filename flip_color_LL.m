function X = flip_color_LL( ...
    X , flips , prior_ll , cell_consts , STA_W , coneConv , colorDot , sizes , beta )
% X = flip_color_LL( X , flips , prior_ll , cell_consts , ...
%                    STA_W , coneConv , colorDot , sizes , beta )
% Apply flips to configuration X, and update log-likelihood of X.
% Bits in X.state are flipped, and the inverse X.invWW of WW is updated.
% The matrix inverse is not recalculated each time:
% block matrix inverse update formulas are used to update X.invWW 
% incrementally, for speed.

% default inverse temperature is 1
if nargin<9
    beta = 1 ;
end

% center of cone RF convolution
Nconv   = ceil(size(coneConv,1)/2) ;

% number of colors
Ncolors = size(colorDot,1) ;

% number of pixels without color
NBW     = length(X.state)/Ncolors ;
if NBW ~= floor(NBW)
    error('number of pixels must be multiple of number of colors!')
end

% initialize inverse if necessary
if ~isfield(X,'invWW') || ~isfield(X,'overlaps')
    temp_X.state     = sparse([],[],[],1,size(STA_W,2),200) ; %false(1,size(STA_W,2)) ;
    temp_X.invWW     = [] ;
    temp_X.overlaps  = [] ;
    temp_X.WW        = [] ;
    temp_X.positions = [] ;
    temp_X.colors    = [] ;
    temp_X.sumLconst = sum(log(cell_consts)) ;
    X = flip_color_LL( temp_X , find(X.state) , prior_ll , cell_consts , ...
                      STA_W , coneConv , colorDot , sizes , beta ) ;
end


% has_delete = 0 ;
% Xold = X ;

% apply all the bit flips to X
for i=1:length(flips)
    k = flips(i) ;
    n = size(X.invWW,1) ;
    
    j = sum(find(X.state)<=k) ;
        
    % block matrix inverse update
    if X.state(k)     % update inverse by deleting jth row/column
        inds = [1:j-1 j+1:n] ;
        X.overlaps = X.overlaps(inds,inds) ;
        X.invWW = X.invWW(inds,inds) - ...
                  X.invWW(inds,j   )*X.invWW(j,inds)/X.invWW(j,j) ;
        X.WW    = X.WW(inds,inds) ;
        X.positions = X.positions(:,inds) ;
        X.colors    = X.colors(inds) ;
        X.state(k)  = false ;
        
%         has_delete = 1 ;
        
        
    else                % update inverse by adding row/column
        j = j + 1 ;
        inds = [1:j-1 j+1:n+1] ;
        A = X.invWW ;
                
%         [kI,kJ] = ind2sub(sizes,mod(k     -1,NBW)+1) ;
        k_position = 1 + mod( k-1 , NBW      ) ;
        kI         = 1 + mod( k-1 , sizes(1) ) ;
        kJ         = 1 + floor( (k_position-1) / sizes(1)) ;
        
        k_color    = 1 + floor( (k-1)          / NBW     ) ;

        % reconstruct W(k,X.state) using coneConv
        Wkk     = coneConv(Nconv,Nconv) ;

        positions   = X.positions ;
        X.positions = zeros(2,n + 1) ;
        if ~isempty(positions)
            X.positions(:,inds) = positions ;
            Wkinds  = [positions(1,:)-kI ; positions(2,:)-kJ] ;
        else
            Wkinds  = [] ;
        end
        X.positions(:,j)    = [kI;kJ] ;
        
        colors      = X.colors ;
        X.colors    = zeros(1,n + 1) ;
        X.colors(inds) = colors ;
        X.colors(j)    = k_color  ;

        overlap = zeros(1,n) ;
        Wkstate = overlap ;
        where   = max(abs(Wkinds),[],1)<Nconv ;
        if sum(where)
            Nc   = size(coneConv,1) ;
            here = (Nconv+Wkinds(2,where) - 1)*Nc + Nconv+Wkinds(1,where) ;
            overlap(where) = coneConv(here) ;
            Wkstate(where) = overlap(where) .* colorDot(k_color,colors(where)) ;
        end

        Wkkc = Wkk * colorDot(k_color,k_color) ;
        
        WW                  = X.WW ;
        X.WW                = zeros(n+1) ;
        X.WW(inds,inds)     = WW ;
        X.WW(inds,j)        = Wkstate ;
        X.WW(j,inds)        = Wkstate ;
        X.WW(j,j)           = Wkkc ;

        O                     = X.overlaps ;
        X.overlaps            = zeros(n+1) ;
        X.overlaps(inds,inds) = O ;
        X.overlaps(inds,j)    = overlap ;
        X.overlaps(j,inds)    = overlap ;
        X.overlaps(j,j)       = Wkk ;
        
        r = Wkstate * A ;
        q = 1/( Wkkc - r*Wkstate' ) ;
        c = A * Wkstate' * q ;
        X.invWW = zeros(n+1) ;
        X.invWW(inds,inds) = A+c*r ;
        X.invWW(inds,j)    = -c    ;
        X.invWW(j,inds)    = -r*q  ;
        X.invWW(j,j)       = q     ;
        
        X.state(k)  = true ;
        
       
        
        
        
%         % debugging!  reconstruct positions from scratch, compare w/
%         % incremental version stored in X
%         this    = find( X.state ) ;
%         posi    = 1 + mod( this-1 , NBW      ) ;
%         tI      = 1 + mod( this-1 , sizes(1) ) ;
%         tJ      = 1 + floor( (posi-1) / sizes(1)) ;        
%         t_color = 1 + floor( (this-1) / NBW     ) ;
%         if      sum(abs(tI - X.positions(1,:)))>0
%             fprintf('case 1')
%         elseif sum(abs(tJ - X.positions(2,:)))>0
%             fprintf('case 2')
%         elseif sum(abs(t_color - X.colors))>0
%             fprintf('case 3')
%         end
        
%         % debugging!  reconstruct WW from scratch
%         overlap = zeros(n+1) ;
%         WW      = zeros(n+1) ;
%         for i=1:n+1
%             for j=1:n+1
%                 dI = min(Nconv-1 , abs(X.positions(1,i) - X.positions(1,j))) + Nconv ;
%                 dJ = min(Nconv-1 , abs(X.positions(2,i) - X.positions(2,j))) + Nconv ;
%                 overlap(i,j) = coneConv( dI , dJ ) ;
%                 WW(i,j)      = coneConv( dI , dJ ) * colorDot(X.colors(j),X.colors(i)) ;
%             end
%         end
%         if      sum(abs(overlap(:) - X.overlaps(:)))>0
%             fprintf('case 4\t')
%             overlap - X.overlaps
%             fprintf('\n')
%         end
%         if      sum(abs(WW(:) - X.WW(:)))>0
%             fprintf('case 5\t')
%             WW - X.WW
%             fprintf('\n')
%         end
%         
%         % debugging!  reconstruct invWW from X.WW
%         invWW = X.WW^(-1) ;
%         if      sum(abs(invWW(:) - X.invWW(:))/max(abs(invWW(:))))>1e-12
%             fprintf('case 6\t')
%             invWW - X.invWW
%             fprintf('\n')
%         end
    end    
end
X.state = logical(X.state) ;

% recalculate data log-likelihood
Ncones = sum(X.state) ;
if Ncones>0
    
    invWW = X.invWW ;
    invWW(abs(invWW)<abs(invWW(1,1))*1e-17) = 0 ;
    invWW = sparse(invWW) ;
    try
        ldet = 2 * sum(log(diag(chol(invWW))));
    catch
        fprintf('\n\nWARNING: catching numerical instability...')
        X.invWW = X.WW^(-1) ;
        ldet = 2 * sum(log(diag(chol(X.invWW))));        
    end
    
    STA_W_state = STA_W(:,X.state) ;

    X.data_ll = full(+ Ncones * (length(cell_consts) * log(2*pi) + X.sumLconst) * ldet + ...
        sum( cell_consts .* sum( (STA_W_state * invWW) .* STA_W_state ,2) )/2) ;
%     X.data_ll = ( - length(cell_consts) * log( det(2.*pi.*X.invWW) ) + ...
%                 sum( cell_consts .* ...
%                       sum( (STA_W(:,X.state) * X.invWW) .* STA_W(:,X.state) ,2) ) ...
%                 )/2 ;
else
    X.data_ll = 0 ;
end

% update log-likelihood
X.ll = beta * (X.data_ll + prior_ll(X)) ;

% if has_delete
%     fprintf(' \t deletion? ,  delta_LL = %.0f' , X.ll - Xold.ll )
% end

% % initialize best 10 configurations encountered to zero
% if ~isfield(X,'best')
%     X.best = zeros(10,size(STA_W,2)+1) ;
%     X.best(1:10) = -Inf ;
% end
% 
% % update best 10 configurations encountered so far
% i = 1 ;
% while 1    
%     if i>10
%         X.best = [ X.best(2:end,:) ; X.ll X.state ] ;
%         break
%     elseif sum(X.state ~= X.best(i,2:end))>0
%         if X.ll > X.best(i)
%             i = i + 1 ;
%         else
%             X.best = [X.best(1:i-1,:) ; X.ll X.state ; X.best(i:10,:)] ;
%             X.best = X.best(2:end,:) ;
%             break
%         end
%     else
%         break
%     end
% end

end