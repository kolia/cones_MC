function [X,done] = propagate_action( X , xy , d , action , done )

% keep track of which cones have already been visited
if nargin<7  , done = sparse([],[],[],X.M0*X.M1,1,10) ;  end

if ~done(xy)
    for cxy = find( X.contact{d}(xy,:) )
        if cxy ~= xy
            [X,done] = propagate_action( X , cxy , d , action , done ) ;
        end
    end
    
    X = action(X,xy) ;
end
done(xy) = true ;

end