function [X,done] = propagate_action( X , x , y , id , d , action , done )

% keep track of which cones have already been visited
if nargin<7  , done = false(X.maxcones,1) ;  end

fprintf('\ndone with %d cones  --  current id %d',sum(done),id)

if ~done(id)
    for cid = find( X.contact{d}(id,:) )
        if cid ~= id
            [X,done] = propagate_action( X , cx , cy , cid , d , action , done ) ;
        end
    end
    
    X = action(X,x,y) ;
end
done(id) = true ;

end