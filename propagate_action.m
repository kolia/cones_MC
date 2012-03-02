function [X,done] = propagate_action( X , xy , d , PROB, T , done )

% keep track of which cones have already been visited
if nargin<6  , done = sparse([],[],[],X.M0*X.M1,1,10) ;  end

if ~done(xy)
    
%     check_X(X)
%     
%     if ~X.state(xy)
%        error('trying to act on nonexistent cell') 
%     end
    
    for cxy = find( X.contact{d}(xy,:) )
        if cxy ~= xy
            [X,done] = propagate_action( X , cxy , d , PROB, T , done ) ;
        end
    end
    
    X = action_LL_shift(X,xy,d,PROB,T) ;
%     X = action(X,xy) ;
end
done(xy) = true ;

end