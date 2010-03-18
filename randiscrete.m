function sample = randiscrete(cumprob,m)
% draw m integers from 1:length(cumprob) from cumulative probability vector
% cumprob, where cumprob must be positive, increasing and cumprob(end) = 1.

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