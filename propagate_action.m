function [X,done] = propagate_action( X , x , y , d , action , done )

% keep track of which cones have already been visited
if nargin<5  , done = false(X.maxcones,1) ;  end

id= X.id(x,y) ;

if ~done(id)
    for cid = find( X.contact{d}(id,:) )
        if cid ~= id
            [X,done] = propagate_action( X , cid , d , action , done ) ;
        end
    end
    
    X = action(X,x,y) ;
end
done(id) = true ;

end