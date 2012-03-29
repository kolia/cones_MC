function svg = plot_LL_ncones( greed , mcmc , cast , cone_map )

if ~isempty(greed)
    [x,y] = get_best( [mcmc ; cast ; {greed}] ) ;
    id = [ones(numel(mcmc),1) ; 2*ones(numel(cast),1) ; 0] ;
else
    [x,y] = get_best( [mcmc ; cast ] ) ;
    id = [ones(numel(mcmc),1) ; 2*ones(numel(cast),1) ] ;
end

y = bits_per_spike( y, cone_map ) ;

minx = min(x) ;
maxx = max(x) ;
inds = (x>=minx) & (x<=maxx) ;

id  = id(inds) ;
y   = y(inds) ;
x   = x(inds) ;

mx = min(x) ;
Mx = max(x) ;

m  = min(y) ;
M  = max(y) ;

id3 = find( y == max(y(id == 1)) ) ;
id( id3 ) = 3 ;

id4 = find( y == max(y(id == 2)), 1 ) ;
id( id4 ) = 4 ;

inds = [1:id3-1 id3+1:id4-1 id4+1:numel(x)-1 id3 id4 numel(x)] ;
x  = x(inds) ;
y  = y(inds) ;
id = id(inds) ;

M1 = max(y(id==1)) ;

mxi = y == max(y) ;

g_ncone = x(end) ;
g_ll = y(end) ;

mcmc_ncone = x(end-2) ;

height = 220 ;

y = y - m ;
y = 200 * (1-y/(M-m)) ;
x = 200 * (x-minx)/(max(x)-minx) ;

svg = sprints('<use xlink:href="#%d" transform="translate(%f %f)"/>\n',id,x,y) ;

ymax = y(find(x==max(x),1)) ;

svg = [sprintf('<line x1="0" x2="%f" y1="0" y2="0" text-anchor="middle" stroke="black" opacity="0.4"/>\n',x(mxi)) ...
       sprintf('<text x="-5" y="0" font-size="12" text-anchor="end" baseline-shift="-45%%">%.5f</text>\n',M) ...
       sprintf('<line x1="%f" x2="%f" y1="%f" y2="%f" text-anchor="middle" stroke="black" opacity="0.4"/>\n',200,200,ymax,height) ...
       sprintf('<text x="%f" y="%f" font-size="12" text-anchor="middle" baseline-shift="-110%%">%d</text>\n',200,height,ceil(Mx)) ...
       sprintf('<line x1="%f" x2="0" y1="%f" y2="%f" text-anchor="middle" stroke="black" opacity="0.4"/>\n',x(end-2),y(end-2),y(end-2)) ...
       sprintf('<text x="-5" y="%f" font-size="12" text-anchor="middle" baseline-shift="-110%%">%.5f</text>\n',y(end-2),M1) ...
       sprintf('<text x="100" y="240" font-size="15" text-anchor="middle" baseline-shift="-100%%" font-weight="bold">     number of cones</text>\n') ...
       sprintf('<g transform="translate(-85 25)rotate(-90)"><text font-size="15" text-anchor="end" font-weight="bold">log posterior (bits/spike)</text></g>\n') ...
       svg] ;

if ~isempty(greed)
    svg = [sprintf('<line x1="0" x2="%f" y1="%f" y2="%f" text-anchor="middle" stroke="black" opacity="0.4"/>\n',x(end),y(end),y(end)) ...
           sprintf('<text x="-5" y="%d" font-size="12" text-anchor="end" baseline-shift="-45%%">%.4f</text>\n',y(end),g_ll) ...
           sprintf('<line x1="%f" x2="%f" y1="%f" y2="%f" text-anchor="middle" stroke="black" opacity="0.4"/>\n',x(end),x(end),y(end),height) ...
           sprintf('<text x="%f" y="%f" font-size="12" text-anchor="middle" baseline-shift="-110%%">%d</text>\n',x(end),height,g_ncone)...
           sprintf('<text x="75" y="-35" font-size="18" text-anchor="middle"><tspan fill="green">Greedy</tspan>, <tspan fill="blue">MCMC</tspan> and <tspan fill="red">CAST</tspan></text>\n') ...
           svg] ;
else
    svg = [sprintf('<line x1="%f" x2="%f" y1="%f" y2="%f" text-anchor="middle" stroke="black" opacity="0.4"/>\n',x(end-2),x(end-2),y(end-2),height) ...
           sprintf('<text x="%f" y="%f" font-size="12" text-anchor="middle" baseline-shift="-110%%">%d</text>\n',x(end-2),height,mcmc_ncone)...
sprintf('<text x="75" y="-35" font-size="18" text-anchor="middle"><tspan fill="blue">MCMC</tspan> and <tspan fill="red">CAST</tspan></text>\n') ...
           svg] ;    
end
   
svg = insert_string(svg,'plot_LL_ncones_stub.svg',-40) ;

if ~isempty(greed)
    save_svg_plot(svg,'LL_ncones.svg')
else
    save_svg_plot(svg,'LL_ncones_nogreed.svg')
end

end