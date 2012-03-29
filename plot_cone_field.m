function svg = plot_cone_field( cc , PROB )

alpha = sqrt(sum(cc.^2,3)) ;
[x,y,alpha] = find( alpha ) ;

[alpha,order] = sort( alpha, 1, 'ascend' ) ;
x = x(order) ;
y = y(order) ;

[M0,M1,~] = size(cc) ;

r = floor( 255*cc(x+M0*(y-1)) ./ alpha ) ;
g = floor( 255*cc(  M0*M1 + x+M0*(y-1)) ./ alpha ) ;
b = floor( 255*cc(2*M0*M1 + x+M0*(y-1)) ./ alpha ) ;

evidence = print_evidence( PROB.NICE ) ;

svg = sprints('<use xlink:href="#1" transform="translate(%f %f)" stroke="rgb(%d,%d,%d)" opacity="%.2f"/>\n', ...
               y,x,r,g,b,alpha) ;
           
scale  = 3 ;
width  = scale*(size(PROB.NICE,2)+ 17) ;
height = scale*(size(PROB.NICE,1)+ 85) ;

fid = fopen('plot_cones_field_stub.svg') ;
svg = sprintf(fread(fid,'*char'),width,height,5,20,width,height,scale,scale,evidence,svg) ;

save_svg_plot(svg,'cones_field.svg')

end