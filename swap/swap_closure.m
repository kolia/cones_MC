function [slowX,fastX,class] = swap_closure( slowX , slowT, fastX , fastT, PROB, class )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations. For this, groups of cells that must be swapped together
% are calculated.
    
N        = numel(slowX.state) ;
Oclass   = class(class<=N) ;
Xclass   = class(class >N) - N ;
M0       = PROB.M0 * PROB.SS ;
class    = [Oclass(:)+N ; Xclass(:)] ;

%     fprintf('\t %d,%d',numel(Xclass),numel(Oclass)) ;

if slowX.maxcones >= slowX.N_cones + numel(Oclass)    && ...
   slowX.maxcones >= fastX.N_cones + numel(Xclass)

        xcx      = 1+mod(Xclass-1,M0) ;
        xcy      = 1+floor((Xclass-1)/M0) ;
        xcc      = slowX.state(Xclass) ;

        ocx      = 1+mod(Oclass-1,M0) ;
        ocy      = 1+floor((Oclass-1)/M0) ;
        occ      = fastX.state(Oclass) ;

        for k=1:length(Xclass)
            slowX = flip_LL( slowX , [xcx(k) xcy(k) 0] , PROB, slowT  ) ;
        end
        for k=1:length(Oclass)
            fastX = flip_LL( fastX , [ocx(k) ocy(k) 0] , PROB, fastT ) ;
        end
        for k=1:length(Oclass)
            slowX = flip_LL( slowX , [ocx(k) ocy(k) occ(k)] , PROB, slowT  ) ;
        end
        for k=1:length(Xclass)
            fastX = flip_LL( fastX , [xcx(k) xcy(k) xcc(k)] , PROB, fastT ) ;
        end
else
    slowX = [] ;
    fastX = [] ;
end

% if N==0
%     X
%     otherX
%     keep
%     numel(classes)
%     classes{1}
%     classes
%     save(sprintf('crash_dump_%d.mat',randi(100)),'X','otherX','keep','classes','PROB')
% end

end