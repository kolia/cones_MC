function X = greedy( X , PROB , update_X )

M0 = PROB.M0 * PROB.SS ;
M1 = PROB.M1 * PROB.SS ;

old_ll  = X.ll ;

if ~isfield(X,'greedy_ll')
    X.greedy_ll = cell(PROB.N_colors,1) ;
%     for c=1:PROB.N_colors
%         X.greedy_ll{c} = old_ll * sparse([],[],[],M0,M1) ;
%     end
end

if ~isfield(X,'changed_x')
    for c=1:PROB.N_colors
        X.greedy_ll{c} = PROB.LL(:,:,c) ;
    end
else
    recalculated = 0 ;
    for i=1:numel(X.changed_x)
        x = X.changed_x(i) ;
        y = X.changed_y(i) ;
        if (x-X.last_x)^2 + (y-X.last_y)^2 > X.D^2
            % propose addition of new cone of each color
            for c=1:PROB.N_colors
                if X.state(x,y) ~= c
                    sample = change_cone( X , [x y c] , PROB , [1 1]) ;
%                     tester = change_cone( X , [X.last_x X.last_y 0] , PROB , [1 1]) ;
%                     tester = change_cone( tester , [x y c] , PROB , [1 1]) ;
                    X.greedy_ll{c}(x,y) = sample.ll - (X.ll + 0.5 * PROB.N_cones_term) ;
%                     if c == 1
%                         fprintf('   %d,%d  %f,%f',x,y,PROB.LL(x,y,c),X.greedy_ll{c}(x,y))
%                         disp('\n')
% %                         fprintf('%d,%d  %f\n',x,y,X.greedy_ll{c}(x,y))
%                     end
                    recalculated = recalculated + 1 ;
                end
            end
        else
            for c=1:PROB.N_colors
                X.greedy_ll{c}(x,y) = -Inf ;
            end
        end
    end
    
%     figure(3)
%     ll = X.greedy_ll{1}(min(X.changed_x)+2:max(X.changed_x)-2,min(X.changed_y)+2:max(X.changed_y)-2) ;
%     imagesc(ll/max(ll(:)))
%     fprintf('changed greedy_ll min %f median %f max %f',min(ll(:)), median(ll(:)), max(ll(:)))
end

% LL = zeros([size(X.greedy_ll{1}) PROB.N_colors]) ;
% for c=1:PROB.N_colors
%     LL(:,:,c) = X.greedy_ll{c} ;
% end
% LL(LL<0) = min(reshape(LL(LL>0),1,[])) ;
% fprintf('   LL min %f max %f',min(LL(:)), max(LL(:)))
% figure(1)
% imagesc(LL/max(LL(:)))
% 
% try
%     figure(2)
%     lll = PROB.LL(min(X.changed_x)+2:max(X.changed_x)-2,min(X.changed_y)+2:max(X.changed_y)-2,1) ;
%     lll( ll == -Inf ) = -Inf ;
%     imagesc(lll/max(lll(:)))
%     fprintf('LL in changed area min %f median %f max %f',min(lll(:)), median(lll(:)), max(lll(:)))
%     disp('done')
% end

% figure(2)
% plot_cones(X.state,PROB) ;

m = -Inf * ones(PROB.N_colors,1) ;
for c=1:PROB.N_colors
    m(c) = max(X.greedy_ll{c}(:)) ;
end

mm = max(m) ;
mc = find(m == mm) ;

if mm>old_ll
    [mx,my] = find(X.greedy_ll{mc} == mm) ;
    
    mx = mx(1) ;
    my = my(1) ;
    
    my = 1+mod(my-1,M1) ;
    
    sx = mod(mx-1,PROB.SS)+1 ;
    sy = mod(my-1,PROB.SS)+1 ;
    [changed_x, changed_y] = find( squeeze(PROB.coneConv(3:end-2,3:end-2,sx,sy)) > 0 ) ;
    
    changed_x = changed_x + mx - sx - PROB.R ;
    changed_y = changed_y + my - sy - PROB.R ;
    
    keep = logical( (changed_x>0) .* (changed_x<=M0) .* (changed_y>0) .* (changed_y<=M1) ) ;
    
    X.changed_x = changed_x(keep) ;
    X.changed_y = changed_y(keep) ;
    X.last_x    = mx ;
    X.last_y    = my ;
    X.last_c    = mc ;
    
    X = change_cone( X , [mx my mc] , PROB , [1 1]) ;
    X = update_X({X}) ;
else
    X = rmfield(X,{'changed_x','changed_y','last_x','last_y'}) ;
end

end