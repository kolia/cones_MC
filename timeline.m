function timeline()

% ( greed , mcmc , cast , cone_map )

n = 1828 ;
gr = [1:n ; cumsum(sort(0.6*rand(1,n),'descend')) ; zeros(1,n) ; zeros(1,n)] ;

mc = cell(30,1) ;
ca = cell(30,1) ;

for i=1:30
    n = ceil(1700+200*rand()) ;
    mc{i} = [1:n ; cumsum(sort((0.9+0.2*rand())*rand(1,n),'descend')) ; zeros(1,n) ; zeros(1,n)] ;
    n = ceil(1700+200*rand()) ;
    ca{i} = [1:n ; cumsum(sort((1.2+0.2*rand())*rand(1,n),'descend')) ; zeros(1,n) ; zeros(1,n)] ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width  = 200 ;
height = 100 ;

xmax = max(gr(1,:)) ;
ymax = max(gr(2,:)) ;
for i=1:30
    xmax = max([xmax mc{i}(1,:)]) ;
    xmax = max([xmax ca{i}(1,:)]) ;
    ymax = max([ymax mc{i}(2,:)]) ;
    ymax = max([ymax ca{i}(2,:)]) ;
end

gr(3:4,:) = [width*gr(1,:)/xmax ; height*gr(2,:)/ymax] ;
for i=1:30
    mc{i}(3:4,:) = [width*mc{i}(1,:)/xmax ; height*mc{i}(2,:)/ymax] ;
    ca{i}(3:4,:) = [width*ca{i}(1,:)/xmax ; height*ca{i}(2,:)/ymax] ;
end

svg = plot_line(gr(3,:),gr(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke="black"/>\n') ;

for i=1:30
    svg = [svg plot_line(mc{i}(3,:),mc{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke="blue"/>\n')] ;
    svg = [svg plot_line(ca{i}(3,:),ca{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke="red"/>\n')] ;
end

fid = fopen('timeline_stub.svg') ;
svg = sprintf(fread(fid,'*char'),svg) ;

save_svg_plot(svg,sprintf('timeline_%s.svg','test'))

end

function svg = plot_line(x,y,pattern,points)
if nargin<4
    points = 20 ;
end
inds = 1:ceil(size(x,2)/points):size(x,2) ;
svg = [] ;
for j=1:length(inds)-1
% svg = [svg sprintf('<use xlink:href="#g" transform="translate(%f %f)"/>\n', ...
%               gr(3,j),gr(4,j))] ;
svg = [svg sprintf(pattern, ...
              x(inds(j)),x(inds(j+1)),-y(inds(j)),-y(inds(j+1)))] ;
end
end