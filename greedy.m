function [X,done] = greedy( X , PROB , update_X )

% if X.N_cones == 0
%     profile clear
%     profile on
% elseif mod(X.N_cones,10) == 0
%     p = profile('info');
%     save(sprintf('profdat_%d',X.N_cones),'p')
%     profile clear
% end

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

oldll = X.ll ;
if ~isfield(X,'changed_x')
    X.greedy_ll = PROB.LL ;
else
    used = 0 ;
    inds = zeros(PROB.N_colors*numel(X.changed_x),1) ;
    gree = zeros(PROB.N_colors*numel(X.changed_x),1) ;
    for i=1:numel(X.changed_x)
        x = X.changed_x(i) ;
        y = X.changed_y(i) ;
        ne = not_excluded( X, x, y ) ;
        for c=1:PROB.N_colors
            used = used + 1 ;
            inds(used) = x + (y-1)*M0 + (c-1)*M0*M1 ;
            if ne && ~isempty(PROB.sparse_struct{x,y,c})
                sample = flip_LL( X , [x y c] , PROB , [1 1]) ;
                gree(used) = sample.ll - X.ll ;
            else
                gree(used) = -Inf ;
            end
        end
    end

    X.greedy_ll(inds(1:used)) = gree(1:used) ;

%     figure(3)
%     ll = X.greedy_ll{1}(min(X.changed_x)+2:max(X.changed_x)-2,min(X.changed_y)+2:max(X.changed_y)-2) ;
%     imagesc(ll/max(ll(:)))
%     fprintf('changed greedy_ll min %f median %f max %f',min(ll(:)), median(ll(:)), max(ll(:)))
end

% if mod( X.N_cones, 20 ) == 1
%     fprintf('   LL min %f max %f',min(X.greedy_ll(isfinite(X.greedy_ll))), max(X.greedy_ll(:)))
%     figure(1)
%     pe = plotable_evidence(X.greedy_ll) ; 
%     imagesc(pe)
% 
%     figure(2)
%     plot_cones(X.state,PROB) ;
%     'word'
% end

[mm,I] = max(X.greedy_ll(:)) ;
[mx,my,mc] = ind2sub(size(PROB.LL),I) ;

done = true ;
if mm>0
    
    newX = flip_LL( X , [mx my mc] , PROB , [1 1]) ;
    if newX.ll>=X.ll
        X = update_X({newX},1,false) ;
        done = false ;
    end

    try
        fprintf('   #keep_cones %d, #keep_GCs %d    mm-dll %f   mm-PROB.ll %f',...
            nnz(X.keep_cones),numel(X.keep_GCs),mm-newX.ll+oldll,mm-PROB.LL(mx,my,mc)) ;
    end


    sx = mod(mx-1,PROB.SS)+1 ;
    sy = mod(my-1,PROB.SS)+1 ;

%     [SIZEX,SIZEY] = size(squeeze(PROB.coneConv(:,:,sx,sy))) ;
%     [changed_x, changed_y] = find( ones(SIZEX+40,SIZEY+40) ) ;
%     changed_x = changed_x-20 ;
%     changed_y = changed_y-20 ;
    [changed_x, changed_y] = find( squeeze(PROB.coneConv(:,:,sx,sy)) > 0 ) ;
    
    changed_x = changed_x + mx - sx - PROB.R ;
    changed_y = changed_y + my - sy - PROB.R ;
    
    keep = logical( (changed_x>0) .* (changed_x<=M0) .* (changed_y>0) .* (changed_y<=M1) ) ;
    
    X.changed_x = changed_x(keep) ;
    X.changed_y = changed_y(keep) ;
    X.last_x    = mx ;
    X.last_y    = my ;
    X.last_c    = mc ;
end

if done
    X = rmfield(X,{'changed_x','changed_y','last_x','last_y'}) ;
end

end