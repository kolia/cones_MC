function timeline( greed , mcmc , cast, fast )

maxtime = max(greed.X.cputime) * 1.1 ;

n_fast = fast.N_cones ;

n  =  greed.X.iteration ;
gr = [greed.X.cputime(1:n)' ; greed.X.LL_history(1:n)' ; zeros(1,n) ; zeros(1,n)] ;

mc = cell(numel(mcmc),1) ;
ca = cell(numel(cast),1) ;

for i=1:numel(mcmc)
    cputimes = reshift_cputime(mcmc{i}.X.cputime,fast) ;
    n = mcmc{i}.X.iteration ;
    cutoff = find(cputimes > maxtime,1) ;
    if ~isempty(cutoff)
        n = min( mcmc{i}.X.iteration, cutoff ) ;
    end
    mc{i} = [cputimes(1:n) ; fast.LL_history(1:n_fast)' mcmc{i}.X.LL_history(n_fast+(1:n-n_fast))' ; zeros(1,n) ; zeros(1,n)] ;
end

for i=1:numel(cast)
    n = find(cast{i}.X.cputime > maxtime,1) ;
    cputimes = reshift_cputime(cast{i}.X.cputime,fast) ;
    ca{i} = [cputimes(1:n) ; fast.LL_history(1:n_fast)' cast{i}.X.LL_history(n_fast+(1:n-n_fast))' ; zeros(1,n) ; zeros(1,n)] ;
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

svg = [] ;

for i=1:numel(mc)
    svg = [svg plot_line(mc{i}(3,:),mc{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.1" stroke="blue"/>\n')] ;
end

for i=1:numel(ca)
    svg = [svg plot_line(ca{i}(3,:),ca{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.1" stroke="red"/>\n')] ;
end

svg = [svg plot_line(gr(3,:),gr(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.1" stroke="black"/>\n')] ;

fid = fopen('timeline_stub.svg') ;
svg = sprintf(fread(fid,'*char'),svg) ;

save_svg_plot(svg,sprintf('timeline__%s.svg','test'))

end

function cputimes = reshift_cputime(cputimes,fast)
m = fast.N_cones ;
cputimes = [fast.cputime(1:m)' cputimes(m+1:end)'-cputimes(m+1)+fast.cputime(m)] - fast.cputime(1) ;
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