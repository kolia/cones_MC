function [X,done] = propagate_action( X , id , d , action , done )

% keep track of which cones have already been visited
if nargin<5  , done = false(X.maxcones,1) ;  end

if ~done(id)
    for cid = find( X.contact{d}(id,:) )
        [X,done] = propagate_action( X , cid , d , action , done ) ;
    end
    
    X = action(X,id) ;
end
done(id) = true ;

end