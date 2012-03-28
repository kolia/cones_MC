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

svg = insert_string([evidence svg],'plot_cones_field_stub.svg',-40) ;

width  = 3*size(PROB.NICE,2)+ 50 ;
height = 3*size(PROB.NICE,1)+ 85 ;

svg = [svg(1:5) sprintf('width="%d" height="%d" viewBox="%d %d %d %d" ',...
                         width,height,5,20,width,height) svg(6:end)] ;

save_svg_plot(svg,'cones_field.svg')

end