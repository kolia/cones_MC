function sample = randiscrete(cumprob,m)
% draw m integers from 1:length(p) from discrete probability vector p.
% p must sum to one.

if nargin<2
    m = 1 ;
end

% cumprob=cumsum(p(:));
sample=zeros(1,m);
for j=1:m
    uni=rand();
    sample(j)=find( uni<=cumprob , 1 ) ;
end

% note
% an implementation with sortrows would be much faster for large m

end