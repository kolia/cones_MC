function X = greedy( X , flip_LL )
% X = greedy( X , flip_LL )
% propagate configuration X through a greedy iteration:
%   - calculate log-likelihood flip_LL of every single flip move
%   - apply to X the move with highest log-likelihood increase

N = length(X.state(:)) ;

% initialize bit location history f necessary
if ~isfield(X,'history')
    h = [] ;
else
    h = X.history ;
end

% initialize output aux
Y0 = flip_LL( X , [] ) ;

% initialize trial configurations y and log-likelihoods
ll      = zeros(1,N+1) ;
ll(1)   = Y0.ll ;
clear X

% calculate the log-likelihoods for each bit flip, and populate y
for i=2:N+1
    Y = flip_LL( Y0 , i-1 ) ;
    ll(i) = Y.ll ;
end

% choose next state
[dummy,i] = max(ll) ;
if i(1)>1
    X = flip_LL( Y0 , i(1)-1 ) ;
else
    X = Y0 ;
end

% imagesc(reshape(X.state,26,46,3))
% drawnow

if ismember(i,h)
    X.history = h( h ~= i-1 ) ;
elseif i>1
    X.history = [h i-1] ;
end

end