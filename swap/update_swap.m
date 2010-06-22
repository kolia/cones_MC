function [results,X] = update_swap(results,X,changed)

% check_X(XLL,X.X)
% check_X(OLL,X.with)

% accumulate acceptance rate statistics
if isfield(results,'N500')
    results.N500     = (1-1/500)*results.N500     + 1 ;
    results.accepted = (1-1/500)*results.accepted + changed ;
%     fprintf('\tX.i=%d,%d:%d,%d:%d',X.X.i,X.with.i,floor(results.N500),floor(results.accepted),changed)
end

[dummy,X.X   ]  = update_X(struct,X.X   ,changed) ;
[dummy,X.with]  = update_X(struct,X.with,changed) ;

% X.X    = xX ;
% X.with = xwith ;

% fprintf('\t swapped ')

% check_X(XLL,X.X)
% check_X(OLL,X.with)

end