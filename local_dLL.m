function dLL = local_dLL( LL , x , y , d , X ) 

c = X.state(x,y) ;

shift = X.masks.shift{d} ;

nx = x + shift(1) ;
ny = y + shift(2) ;

if nx>0 && nx<=X.M0 && ny>0 && ny<=X.M1
    dLL = LL(nx,ny,c)-LL(x,y,c) ;
else
    dLL = -LL(x,y,c) + X.N_cones_factor ;
end

end