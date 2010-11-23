function svg = plot_G_MC_PT( greedy , dir1 , pattern1 , dir2 , pattern2 )

phase1 = load_bestX( dir1 , pattern1 ) ;
phase2 = load_bestX( dir2 , pattern2 ) ;

ll = get_ll( [phase1 ; phase2 ; {greedy}] ) ;
id = [ones(numel(phase1),1) ; 2*ones(numel(phase2),1) ; 0] ;

[xG,yG,cG] = find( greedy.state ) ;
[x1,y1,c1] = find( phase1{ ll == max(ll(id == 1)) }.state ) ;
[x2,y2,c2] = find( phase2{ find( ll == max(ll(id == 2)) , 1) - numel(phase1) }.state ) ;

id = [ones(numel(x1),1) ; 2*ones(numel(x2),1) ; zeros(numel(xG),1)] ;

c  = [num2cell(c1(:)) ; num2cell(c2(:)) ; num2cell(cG(:))] ;
for i=1:numel(c)
    cc = c{i} ;
    switch cc
        case 1
            c{i} = 'red' ;
        case 2
            c{i} = 'green' ;
        case 3
            c{i} = 'blue' ;
    end
end

svg = sprints('<use xlink:href="#%d" transform="translate(%f %f)" stroke="%s"/>\n',id,[y1(:);y2(:);yG],[x1(:);x2(:);xG],c) ;

% svg = [sprintf('<text x="100" y="220" font-size="15" text-anchor="middle" baseline-shift="-100%%">     number of cones</text>\n') ...
%        sprintf('<g transform="translate(-10 55)rotate(-90)"><text font-size="15" text-anchor="end">log posterior</text></g>\n') ...
%        sprintf('<text x="100" y="-25" font-size="18" text-anchor="middle"><tspan fill="green">Greedy</tspan>, <tspan fill="blue">MCMC</tspan> and <tspan fill="red">Parallel tempering</tspan></text>\n') ...
%        svg] ;

svg = insert_string(svg,'plot_G_MC_PT_stub.svg',-40) ;

fid = fopen('plot_G_MC_PT.svg','w') ;
fwrite(fid,svg) ; fclose(fid) ;

end


function y = get_ll( Xs )

y = zeros(numel(Xs),1) ;
for i=1:numel(Xs)
    y(i) = Xs{i}.ll ;
end

end