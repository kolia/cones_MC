function svg = plot_Greedy_MCMC_CAST( greed , mcmc , cast , NICE )

[~,ll,states,keep] = get_best( [mcmc ; cast ; {greed}] ) ;
id = [ones(numel(mcmc),1) ; 2*ones(numel(cast),1) ; 0] ;
id = id(keep) ;

[xG,yG,cG] = find( states{end} ) ;
[x1,y1,c1] = find( states{ ll(id == 1) == max(ll(id == 1)) } ) ;
[x2,y2,c2] = find( states{ ll(id == 2) == max(ll(id == 2)) } ) ;

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

evidence = print_evidence( NICE ) ;

svg = sprints('<use xlink:href="#%d" transform="translate(%f %f)" stroke="%s"/>\n', ...
               id,[y1(:);y2(:);yG],[x1(:);x2(:);xG],c) ;

svg = insert_string([evidence svg],'plot_Greedy_MCMC_CAST_stub.svg',-40) ;
save_svg_plot(svg,'Best_Greed_MCMC_CAST.svg')

end