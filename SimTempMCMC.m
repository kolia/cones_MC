function [X,ST] = SimTempMCMC( X, PROB, ST , j )
% Simulated Tempering + Wang-Landau update
% simTempMCMC : We update the state, and then update the temperature
%   ST (Simulated Tempering) is a structure that contains fields:
%       T  : cell array of temperatures
%       i  : index of current temperature
%       lg : current log weights
%       n  : number of current iteration
    
if ~isfield(ST,'n'     ) , ST.n      = 1 ;                       end
if ~isfield(ST,'i'     ) , ST.i      = 1 ;                       end
if ~isfield(ST,'lg'    ) , ST.lg     = zeros(length(ST.T),1)   ; end
if ~isfield(ST,'f'     ) , ST.f      = zeros(length(ST.T),1)   ; end
if ~isfield(ST,'gamma0') , ST.gamma0 = 1e128 ; end
if ~isfield(ST,'k'     ) , ST.k      = 0 ; end
if ~isfield(ST,'avg_p' ) , ST.avg_p  = zeros(length(ST.T),1) ; end
if ~isfield(ST,'Q'     )
    ST.Q = ones(length(ST.T),1) ;
    ST.Q(2:end-1) = 0.5 ;
end

ST.gamma   = ST.gamma0 ;
ST.lg_curv = 0.3 * ones(length(ST.T)-1,1) ;

% sample T
% proposed change in temperature index is +1 or -1 with prob 0.5
if ST.i(j) == 1
    proposed_i = 2 ;
elseif ST.i(j) == length(ST.T)
    proposed_i = length(ST.T)-1 ;
else
    proposed_i = ST.i(j) + 2* ( unifrnd(0,1)>0.5 ) - 1 ;
end

% calculate probability of accepting new T
new_ll  = calculate_LL( X , PROB , ST.T{proposed_i}) ;
paccept = min(1, ST.Q(ST.i(j))/ST.Q(proposed_i)*...
                 exp(new_ll-X.ll + ST.lg(ST.i(j)) - ST.lg(proposed_i)) ) ;

ST.avg_p( ST.i(j) ) = ST.avg_p( ST.i(j) ) + paccept ;

% accept new state with probability p
accept = rand() < paccept ;
if accept
    ST.i(j) = proposed_i ;
    X.ll = new_ll ;
end

ST.f(ST.i(j)) = ST.f(ST.i(j)) + 1 ;

% update adaptive weights
if numel(ST.f)*ST.f/sum(ST.f)>0.8
%     fprintf('\nST.lg:') ; fprintf(' %g',ST.lg)
%     fprintf('\nST.f:') ;  fprintf(' %g',ST.f/s)
    ST.k = ST.k+1 ;
    ST.gamma = (1 + ST.gamma0).^(1/(2*ST.k)) - 1 ;
%     fprintf('\navg_acceptance') ; fprintf(' %g',ST.avg_p./ST.f) ;
    
    if isfield(ST,'max_temps')  &&  2*length(ST.T)-1 <= ST.max_temps
        [ST.T,T] = refine_mesh(ST.T    , ST.curvature) ;
        ST.lg    = spline(1:length(ST.lg),ST.lg,1:0.5:length(ST.lg)) ;
        ST.f     = zeros(length(ST.T),1) ;
        ST.avg_p = zeros(length(ST.T),1) ;
        ST.Q     = ones( length(ST.T),1) ;
        
        ST.curvature = 1./(1+sqrt((1-ST.curvature)/ST.curvature)) ;

        ST.Q(2:end-1)= 0.5 ;
        ST.i = 2*ST.i-1 ;
        
        fprintf('\n\n+ %d inverse temperatures beta:\n',size(T,1))
        fprintf('%4.2f  ',T(:,1)' )
        fprintf('\n\n+ %d powers delta:\n',size(T,1))
        fprintf('%4.2f  ',T(:,2)')
    end
    fprintf('\nNEW GAMMA: %g  (k=%d)\n',ST.gamma,ST.k)
end

ST.lg(ST.i(j)) = ST.lg(ST.i(j)) + log(1 + ST.gamma) ;

ST.n       = ST.n + 1 ;

% reshift ST.lg every 1000 iterations, for sanity
if ~mod(ST.n,1000) ,  ST.lg = ST.lg - mean(ST.lg) ; end

end