function [results,X] = update_swap(results,X,changed)

% check_X(XLL,X.X)
% check_X(OLL,X.with)

% accumulate acceptance rate statistics
if isfield(results,'N500')
%     fprintf('\tX.i=%d,%d',X.X.i,X.with.i)
    results.N500     = (1-1/500)*results.N500     + 1 ;
    results.accepted = (1-1/500)*results.accepted + changed ;
end

[dummy,xX   ]  = update_X(struct,X.X   ,changed) ;
[dummy,xwith]  = update_X(struct,X.with,changed) ;

X.X    = xX ;
X.with = xwith ;

% fprintf('\t swapped ')

% check_X(XLL,X.X)
% check_X(OLL,X.with)

end