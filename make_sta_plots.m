function make_sta_plots( sta )
% Given a cell struct of stas, for example as output by denoised_sta, save
% imagesc plots to current directory.  Negative regions in STA are overlaid
% with red dots.

NGC = numel(sta) ;
[N,M] = size( sta{1} ) ;

chec = 1:20:N*M ;
checker = zeros(N,M) ;
checker(chec) = 1 ;

for i=1:NGC
    st = sta{i} ;
    z = sum(st,3) ;
    z = z<0 ;
    z = z & checker ;
    st = st-min(st(:)) ;
    st = st / max(st(:)) ;
    st( z ) = 1 ;
    imagesc(st) ;
    saveas(h,sprintf('cell%d',i),'jpg') ;
end

end