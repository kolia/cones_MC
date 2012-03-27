function svg = print_evidence( NICE )

[nx,ny,~] = size(NICE) ;
f = figure('visible','off'); imshow(NICE, 'Border', 'tight');
set(gca,'position',[0 0 1 1],'units','normalized')
set(gcf,'PaperPosition',[0 0 ny/80 nx/80])
set(gcf,'PaperSize',[ny nx])
print(f, '-r80', '-dpng', 'evidence80.png');

svg = sprintf('<image width="%d" height="%d" xlink:href="evidence80.png"/>\n',ny,nx) ;

end