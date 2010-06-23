function [results,X] = update_swap(results,trials,i)

% check_X(XLL,X.X)
% check_X(OLL,X.with)

% accumulate acceptance rate statistics
if isfield(results,'N500')
    results.N500     = (1-1/500)*results.N500     + 1 ;
    results.accepted = (1-1/500)*results.accepted + (i>1) ;
%     fprintf('\tX.i=%d,%d:%d,%d:%d',X.X.i,X.with.i,floor(results.N500),floor(results.accepted),changed)
end

X.X         = trials{1}.X ;
X.X.state   = trials{i}.X.state ;
X.X.invWW   = trials{i}.X.invWW ;
X.X.N_cones = trials{i}.X.N_cones ;
X.X.ll      = trials{i}.X.ll ;
X.X.diff    = trials{i}.X.diff ;
X.X.beta    = trials{i}.X.beta ;

X.with         = trials{1}.with ;
X.with.state   = trials{i}.with.state ;
X.with.invWW   = trials{i}.with.invWW ;
X.with.N_cones = trials{i}.with.N_cones ;
X.with.ll      = trials{i}.with.ll ;
X.with.diff    = trials{i}.with.diff ;
X.with.beta    = trials{i}.with.beta ;

[dummy,X.X   ]  = update_X(struct,{X.X   },1) ;
[dummy,X.with]  = update_X(struct,{X.with},1) ;

if i>1
    X.X.version     = X.X.version    + 1 ; 
    X.with.version  = X.with.version + 1 ;
end

results.trials = trials ;
results.version = [X.X.version X.with.version] ;

% X.X    = xX ;
% X.with = xwith ;

% fprintf('\t swapped ')

% check_X(XLL,X.X)
% check_X(OLL,X.with)

end