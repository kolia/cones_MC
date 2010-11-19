function svg = plot_LL_ncones( greedy , dir1 , pattern1 , dir2 , pattern2 )

phase1 = load_bestX( dir1 , pattern1 ) ;
phase2 = load_bestX( dir2 , pattern2 ) ;

[x,y,m,M] = get_N_ll( [{greedy} ; phase1 ; phase2] ) ;
id = [0 ; ones(numel(phase1),1) ; 2*ones(numel(phase2),1)] ;

svg = svg_use( x , y , id ) ;
svg = [sprintf('<text id="min" x="-10" y="100" >0</text>\n') ...
       sprintf('<text id="min" x="-30" y="0" >%d</text>\n',ceil(M-m)) ...
       svg] ;
   
svg = insert_string(svg,'plot_LL_ncones_stub.svg',-40) ;

fid = fopen('LL_ncones.svg','w') ;
fwrite(fid,svg) ; fclose(fid) ;

end


function [x,y,ymin,ymax] = get_N_ll( Xs )

x = zeros(numel(Xs),1) ;
y = zeros(numel(Xs),1) ;
for i=1:numel(Xs)
    x(i) = Xs{i}.N_cones ;
    y(i) = Xs{i}.ll ;
end
ymin = min(y) ;
ymax = max(y) ;
y = y - ymin ;
y = 100 * (1-y/max(y)) ;

end