function [X1,X2] = swap_step(X1,T1,X2,T2,PROB)

LL = @(x)get_LL(x.X,PROB,T1)+get_LL(x.with,PROB,T2) ;
% check_X(X1)
% check_X(X2)
swapX = swap_closure( X1, T1, X2 , T2, PROB) ;
% for i=1:numel(swapX)
%     check_X(swapX{i}.X)
%     check_X(swapX{i}.with)
% end
swapX = flip_MCMC( swapX{1}, swapX(2:end), @update_swap , LL ) ;
X1  = swapX.X ;
X2  = swapX.with ;

end
