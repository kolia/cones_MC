function samples = independent_sampler( X , p , n )
% produce n independent samples from discrete distribution p
% in a format legible by flip_MCMC - a cell structure with fields:
%   - samples{i}.flips, # of the flipped bit for sample i
%   - samples{i}.forward_prob , the prob of choosing this trial move
%   - samples{i}.backward_prob, the prob of choosing the reverse move

sampled = randiscrete( p , n) ;
samples = cell(n,1) ;
for i=1:n
    samples{i}.flips         = sampled(i) ;
    samples{i}.forward_prob  = p(sampled(i)) ;
    samples{i}.backward_prob = p(sampled(i)) ;
end

end