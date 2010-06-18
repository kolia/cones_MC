function results = accumulate_results(results,X,i)

results.summed = results.summed + ...
    [(X.state(:)'==1) (X.state(:)'==2) (X.state(:)'==3)] ;

% accumulate acceptance rate statistics
if isfield(results,'stats')
    results.stats.N500 = (1-1/500)*results.stats.N500 + 1 ;
    results.stats.accepted = (1 - 1/500) * results.stats.accepted + (i>1) ;
end

end