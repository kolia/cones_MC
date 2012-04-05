function [free,indices] = not_excluded(X,x,y)

[~,indices] = place_mask( X.M0 , X.M1 , x , y , X.masks{1,1}.exclusion ) ;
free    = isempty( find( X.state(indices)>0 , 1) ) ;

end
