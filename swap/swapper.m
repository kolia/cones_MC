function swap = swapper( X , otherX )
% Combine two systems into a single system, summing their log-likelihoods.
% This is in preparation for swapping parts of the two systems'
% configurations:  connected components of the symmetric difference between
% the two configurations are calculated and stored in swap.flips.

swap.X = X ;
swap.with = otherX ;

sym_diff = abs(X.state - otherX.state) ;
x = zeros(length(X.state),1) ;

% try newer version bwconncomp
try
    CC = bwconncomp( full(sym_diff) , 8 ) ;
   
%     fprintf('\n%d swaps: ')
    
    swap.flips = cell(CC.NumObjects,1) ;
    for j=1:CC.NumObjects
        
%         fprintf('%d ',length(CC.PixelIdxList{j}))
        
        swap.flips{j}.flips = CC.PixelIdxList{j} ;
        swap.flips{j}.forward_prob  = 1/CC.NumObjects ;
        swap.flips{j}.backward_prob = 1/CC.NumObjects ;
        
        x(CC.PixelIdxList{j}) = j ;
    end
       
% or fall back to older bwlabel
catch
    [L, num] = bwlabel(full(sym_diff), 8) ;
    swap.flips = cell(num,1) ;
    for j=1:num
        swap.flips{j}.flips = find(L == j) ;
        swap.flips{j}.forward_prob  = 1/num ;
        swap.flips{j}.backward_prob = 1/num ;
        
        x( swap.flips{j}.flips ) = j ;
    end

end
    
% imagesc(x)
% title('Connected components','FontSize',20)
% colorbar
% drawnow
% waitforbuttonpress

swap.ll = X.ll + otherX.ll ;
swap.state = X.state ;

end