function ST = SimTempMCMC( X, PROB, LL, ST )
%simTempMCMC : We update the state, and then update the temperature
%   ST (Simulated Tempering) is a structure that contains fields:
%       T : cell array of temperatures
%       i : index of current temperature
%       g : current weights
%       n : number of current iteration
    
if ~isfield(ST,'n') , ST.n = 1 ;                    end
if ~isfield(ST,'i') , ST.i = length(ST.T) ;         end
if ~isfield(ST,'g') , ST.g = ones(length(ST.T),1) ; end
    

%step 4 (sample T)
% proposed change in temperature index is +1 or -1 with prob 0.5
di = 2* ( unifrnd(0,1)>0.5 ) - 1 ;
% make sure proposed_i is between 1 and length(ST.T)
proposed_i = min( max( 1 , ST.i + di ) , length(ST.T) ) ;

%calculate probability of accepting new T
paccept = min(1, LL(X, PROB, ST.T{proposed_i})* g(ST.i      ) ...
               /(LL(X, PROB,         T       )* g(proposed_i)));

%accept new state with probability p
if unifrnd(0,1) < paccept
    ST.i = proposed_i ;
end

%step 5 update adaptive weights  MAY NEED TO MAKE A LIST INSTEAD OF A
%NUMBER (use log g)
%logG = log(g);
%logGnew = logG;
%logGnew(ind) = logG(ind) + log(1+gamma);
%gnew(ind) = exp(logGnew(ind));
ST.g(ST.i) = ST.g(ST.i)*(1 + 1/sqrt(1+ST.n));
ST.n       = ST.n + 1 ;

end