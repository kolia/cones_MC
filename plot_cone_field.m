function svg = plot_cone_field( init, cast , PROB )

[x,yc,v] = find( cc ) ;

M1     = size(cc,2) ;
y      = yc(1:M1) ;
colors = floor( 255*reshape(v, [], 3) ) ;

evidence = print_evidence( PROB.NICE ) ;

svg = sprints('<use xlink:href="#1" transform="translate(%f %f)" stroke="rgb(%d,%d,%d)"/>\n', ...
               x,y,colors(:,1),colors(:,2),colors(:,3)) ;

svg = insert_string([evidence svg],'plot_cones_field_stub.svg',-40) ;
save_svg_plot(svg,'cones_field.svg')

end