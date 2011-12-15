function check_LL( Y , PROB )

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

[xs , ys , cs] = find(Y.state) ;

posX    = [] ;
posY    = [] ;
colors  = [] ;

X.state = Y.state * 0 ;
X.N_cones = 0 ;
X.WW    = [] ;

order = randperm(size(xs,1)) ;

for i=order
    
    x = xs(i) ;
    y = ys(i) ;
    c = cs(i) ;
    
    % block matrix inverse update
    X.N_cones   = X.N_cones + 1 ;
    
    Wkinds  = [posX'-x ; posY'-y] ;
    ssx     = 1+mod(x-1,PROB.SS) ;
    ssy     = 1+mod(y-1,PROB.SS) ;
            
    Wkstate = zeros(1,X.N_cones-1) ;

%     where   = find( max(abs(Wkinds),[],1) <= PROB.R ) ;
%     if ~isempty(where)
%         xx = Wkinds(1,where)+PROB.R+ssx ;
%         yy = Wkinds(2,where)+PROB.R+ssy ;
%         for k=1:length(where)
%             Wkstate(where(k)) = PROB.coneConv(xx(k),yy(k),ssx,ssy) ...
%                 .* PROB.colorDot(c,colors(where(k))) ;
%         end
        
    where   = find( max(abs(Wkinds),[],1) <= PROB.R ) ;
    if ~isempty(where)
        xx = Wkinds(1,where)+PROB.R+ssx ;
        yy = Wkinds(2,where)+PROB.R+ssy ;
        for k=1:length(where)
            Wkstate(where(k)) = PROB.coneConv(xx(k),yy(k),ssx,ssy) ...
                .* PROB.colorDot(c,colors(where(k))) ;
        end
        
        
%     for k=1:size(posX)
%         Wkstate(k) = 1/norm(Wkinds(:,k)) ;
%     end
    end
    
    Wkkc    = PROB.coneConv(PROB.R+ssx,PROB.R+ssy,ssx,ssy) * PROB.colorDot(c,c) ;
%     Wkkc    = i ;
    
    if max(Wkstate) > Wkkc
        error('This is absurd.')
    end
    
    X.WW            = [X.WW Wkstate' ; Wkstate Wkkc ] ;
    
    X.state(x,y)    = c ;
    posX            = [posX ; x] ;
    posY            = [posY ; y] ;
    colors          = [colors  ; c] ;
    
end

% recalculate data log-likelihood
if X.N_cones>0
    
    STA_W_state = PROB.STA_W( posX'+M0*(posY'-1)+M0*M1*(colors'-1) , : )' ;
    
    ll  = Y.beta * full(- X.N_cones * PROB.sumLconst - log(det(X.WW)) * PROB.N_GC + ...
        sum( PROB.cell_consts .* sum( (STA_W_state * inv(X.WW)) .* STA_W_state ,2) )/2) ;
else
    ll = 0 ;
end

X.ll  = ll ;

if max(abs(x-xs))>0 || max(abs(y-ys))>0 || max(abs(c-cs))>0
    'problem with check_LL...'
end

% if abs(ll-Y.ll)>1e-8
if abs(det(X.WW) - 1/det(Y.invWW))>1e-8
    'det(WW) is not right'
end

order ;
order = full( sparse(order,1:numel(order),1:numel(order),numel(order),numel(order)) ) ;
order = max(full(order')) ;
[ll-Y.ll ll Y.ll]
[posX(order) posY(order) colors(order)]
[X.WW(order,order) ; inv(Y.invWW)]

if abs(ll - Y.ll)>1e-8
    'Y.ll is not right'
end