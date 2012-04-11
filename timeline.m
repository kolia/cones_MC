function timeline( greed , mcmc , cast, fast, total_spikes )

maxtime = max(greed.X.cputime) * 1.15 ;

ll_fast = fast.ll ;
cpustart = fast.cputime(fast.N_cones) ;

n     = greed.X.iteration ;
start = find( greed.X.LL_history>ll_fast, 1 ) ;
gr = [greed.X.cputime(start:n)' ; greed.X.LL_history(start:n)' ; zeros(1,n-start+1) ; zeros(1,n-start+1)] ;

mc = cell(numel(mcmc),1) ;
ca = cell(numel(cast),1) ;

for i=1:numel(mcmc)
    cputimes = reshift_cputime(mcmc{i}.X.cputime,fast) ;
    n = mcmc{i}.X.iteration ;
    cutoff = find(cputimes > maxtime,1) ;
    if ~isempty(cutoff)
        n = min( mcmc{i}.X.iteration, cutoff ) ;
    end
    start = find( mcmc{i}.X.LL_history>ll_fast, 1) ;
    mc{i} = [cputimes(start:n) ; mcmc{i}.X.LL_history(start+(1:n-start+1))' ; zeros(1,n-start+1) ; zeros(1,n-start+1)] ;
end

maxcast = 0 ;
for i=1:numel(cast)
    maxcast = max(maxcast,max(cast{i}.X.LL_history)) ;
    n = find(cast{i}.X.cputime > maxtime,1) ;
    cputimes = reshift_cputime(cast{i}.X.cputime,fast) ;
    start = find( mcmc{i}.X.LL_history>ll_fast, 1) ;
    ca{i} = [cputimes(start:n) ; cast{i}.X.LL_history(start+(1:n-start+1))'  ; zeros(1,n-start+1) ; zeros(1,n-start+1)] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width  = 230 ;
height = 100 ;

xmax = max(gr(1,:)) ;
ymax = max(gr(2,:)) ;
for i=1:30
    xmax = max([xmax mc{i}(1,:)]) ;
    xmax = max([xmax ca{i}(1,:)]) ;
    ymax = max([ymax mc{i}(2,:)]) ;
    ymax = max([ymax ca{i}(2,:)]) ;
end

ymax = maxcast+100 ;

maxcast_rescaled = height*(maxcast - ll_fast)/(ymax-ll_fast) ;
gr(3:4,:) = [width*gr(1,:)/(xmax) ; height*(gr(2,:) - ll_fast)/(ymax-ll_fast)] ;
for i=1:30
    mc{i}(3:4,:) = [width*mc{i}(1,:)/(xmax) ; height*(mc{i}(2,:) - ll_fast)/(ymax-ll_fast)] ;
    ca{i}(3:4,:) = [width*ca{i}(1,:)/(xmax) ; height*(ca{i}(2,:) - ll_fast)/(ymax-ll_fast)] ;
end

svg = [] ;
svg = [svg sprintf('<text x="-12" y="-3" fill="black" font-size="7" baseline-shift="-35%%" font-family="Helvetica">%.4f</text>\n', bits_per_spike( ll_fast, total_spikes))] ;
svg = [svg sprintf('<text x="%f" y="7" text-anchor="start" fill="black" font-size="7" baseline-shift="-35%%" font-family="Helvetica">%d sec.</text>\n', mc{1}(3,1), floor(cpustart))] ;
svg = [svg sprintf('<text x="%f" y="7" text-anchor="middle" fill="black" font-size="7" baseline-shift="-35%%" font-family="Helvetica">%d sec.</text>\n', gr(3,1), floor(gr(1,1)))] ;
svg = [svg sprintf('<text x="%f" y="%f" fill="black" font-size="7" baseline-shift="-35%%" font-family="Helvetica">%.4f</text>\n',mc{1}(3,1),-maxcast_rescaled,bits_per_spike( maxcast, total_spikes ))] ;
svg = [svg sprintf('<line x1="%f" x2="230" y1="%f" y2="%f" stroke-width="0.6" stroke="red" stroke-dasharray="9 5"/>\n',mc{1}(3,1)+25,-maxcast_rescaled,-maxcast_rescaled)] ;
svg = [svg sprintf('<circle r="1" cx="%f" cy="0"  fill="green"/>\n',mc{1}(3,1))] ;
svg = [svg sprintf('<circle r="1" cx="%f" cy="%f" fill="black"/>\n',gr(3,1),-gr(4,1))] ;
svg = [svg sprintf('<circle r="1" cx="%f" cy="%f" fill="black"/>\n',gr(3,end),-gr(4,end))] ;

for i=1:numel(mc)
    svg = [svg plot_line(mc{i}(3,:),mc{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.3" stroke="blue"/>\n')] ;
end

for i=1:numel(ca)
    svg = [svg plot_line(ca{i}(3,:),ca{i}(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.3" stroke="red"/>\n')] ;
end

svg = [svg plot_line(gr(3,:),gr(4,:),'<line x1="%f" x2="%f" y1="%f" y2="%f" stroke-width="0.6" stroke="black"/>\n')] ;

fid = fopen('timeline_stub.svg') ;
svg = sprintf(fread(fid,'*char'), svg) ;

save_svg_plot(svg,sprintf('timeline__%s.svg','test'))

end

function cputimes = reshift_cputime(cputimes,fast)
m = fast.N_cones ;
cputimes = [fast.cputime(1:m)' cputimes(m+1:end)'-cputimes(m+1)+fast.cputime(m)] - fast.cputime(1) ;
end

function svg = plot_line(x,y,pattern,points)
if nargin<4
    points = 30 ;
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