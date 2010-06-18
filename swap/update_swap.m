function X = update_swap(oldX,X)

% check_X(XLL,X.X)
% check_X(OLL,X.with)

X.X     = update_X(oldX,X.X   ) ;

X.with  = update_X(oldX,X.with) ;

% check_X(XLL,X.X)
% check_X(OLL,X.with)

end