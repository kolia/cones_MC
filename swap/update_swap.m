function [results,X] = update_swap(results,trials,i)

% if isfield(results,'results')       && ...
%    isfield(results.results{1},'i')  && results.results{1}.i == 1
%     check_X(trials{1}.X)
%     check_X(trials{1}.with)
% end

% accumulate acceptance rate statistics
if isfield(results,'N500')
    results.N500     = (1-1/500)*results.N500     + 1 ;
    results.accepted = (1-1/500)*results.accepted + (i>1) ;
%     fprintf('\tX.i=%d,%d:%d,%d:%d',X.X.i,X.with.i,...
%               floor(results.N500),floor(results.accepted),changed)
end

if isfield(results,'N50')
    results.N50      = (1-1/50)*results.N50     + 1 ;
    results.accepted = (1-1/50)*results.accepted + (i>1) ;
%     fprintf('\tX.i=%d,%d:%d,%d:%d',X.X.i,X.with.i,...
%               floor(results.N50),floor(results.accepted),changed)
end

if i>1

%     check_X(trials{1}.X)
%     check_X(trials{1}.with)

    X.X         = trials{1}.X ;
    X.X.state   = trials{i}.X.state ;
    X.X.invWW   = trials{i}.X.invWW ;
    X.X.N_cones = trials{i}.X.N_cones ;
    X.X.ll      = trials{i}.X.ll ;
    X.X.diff    = trials{i}.X.diff ;
    X.X.beta    = trials{i}.X.beta ;
    X.X.delta   = trials{i}.X.delta ;
    
    X.with         = trials{1}.with ;
    X.with.state   = trials{i}.with.state ;
    X.with.invWW   = trials{i}.with.invWW ;
    X.with.N_cones = trials{i}.with.N_cones ;
    X.with.ll      = trials{i}.with.ll ;
    X.with.diff    = trials{i}.with.diff ;
    X.with.beta    = trials{i}.with.beta ;
    X.with.delta   = trials{i}.with.delta ;
    
%     check_X(X.X)
%     check_X(X.with)
    
    X.X.version     = X.X.version    + 1 ;
    X.with.version  = X.with.version + 1 ;
else
    X = trials{1} ;
end

% update both X
ii = 2*(i>1) + (i==1) ;
[results.results{1},X.X   ]  = update_X(results.results{1},{trials{1}.X    X.X   },ii) ;
[results.results{2},X.with]  = update_X(results.results{2},{trials{1}.with X.with},ii) ;    

if isfield(results.results{1},'swap')
    results.results{1}.swap(results.results{1}.iteration) = true ;
end

if isfield(results.results{2},'swap')
    results.results{2}.swap(results.results{2}.iteration) = true ;
end

% if i==1
%     results.trials = trials ;
%     results.version = [X.X.version X.with.version] ;
% end

% if isfield(results,'results')       && ...
%    isfield(results.results{1},'i')  && results.results{1}.i == 1
%     check_X(X.X)
%     check_X(X.with)
% end

end