function svg = plot_cone_field( cc , PROB , filename )
% plot_cone_field( X, cone_map )
% plot_cone_field( X, cone_map, filename )
% plot_cone_field( confidence, cone_map )
% plot_cone_field( confidence, cone_map, filename )

if nargin<3
    filename = 'cones_field.svg' ;
end

if isstruct( cc )
    X = cc ;
    cc = zeros(PROB.M0 * PROB.SS, PROB.M1 * PROB.SS, 3) ;
    for c=1:3
        cc(:,:,c) = cc(:,:,c) + (X.state==c) ;
    end
end

alpha = sqrt(sum(cc.^2,3)) ;
[x,y,alpha] = find( alpha ) ;

[alpha,order] = sort( alpha, 1, 'ascend' ) ;
x = x(order) ;
y = y(order) ;

[M0,M1,~] = size(cc) ;

r = floor( 255*cc(x+M0*(y-1)) ./ alpha ) ;
g = floor( 255*cc(  M0*M1 + x+M0*(y-1)) ./ alpha ) ;
b = floor( 255*cc(2*M0*M1 + x+M0*(y-1)) ./ alpha ) ;

svg = sprints('<use xlink:href="#1" transform="translate(%f %f)" stroke="rgb(%d,%d,%d)" opacity="%.2f"/>\n', ...
               y,x,r,g,b,alpha) ;

scale  = 500/max([size(PROB.NICE,1) size(PROB.NICE,2)]) ;
width  = min([500 500*size(PROB.NICE,2)/size(PROB.NICE,1)])+20 ;
height = min([500 500*size(PROB.NICE,1)/size(PROB.NICE,2)])+20 ;

% width  = scale*(size(PROB.NICE,2)+ 17) ;
% height = scale*(size(PROB.NICE,1)+ 85) ;

% scale = 1 ;
% [width,height] = size(PROB.NICE) ;

evidence = print_evidence( PROB.NICE ) ;
% evidence = [sprintf('<g transform="scale(%f,%f)">',scale,scale) evidence '</g>'] ;

fid = fopen('plot_cones_field_stub.svg') ;
svg = sprintf(fread(fid,'*char'),width,height,width,height,scale,scale,evidence,svg) ;

save_svg_plot(svg,filename)

end