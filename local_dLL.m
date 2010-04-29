function dLL = local_dLL( LL , x , y , shift , M0 , M1 )

nx = x + shift(1) ;
ny = y + shift(2) ;

if nx>0 && nx<=M0 && ny>0 && ny<=M1
    dLL = LL(nx,ny)-LL(x,y) ;
else
    dLL = -LL(x,y) ;
end

end