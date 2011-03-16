function make_sta_plots( sta , filenames )
% Given a cell struct of stas, for example as output by denoised_sta, save
% imagesc plots to current directory.  Negative regions in STA are overlaid
% with red dots.

NGC = numel(sta) ;
[N,M] = size( sta{1} ) ;

h = figure ;

for i=1:NGC
    st = sta{i} ;
    z = sum(st,3) ;
    inds  = logical(z<0) ;
    inds  = inds(:) ;
    inds3 = [inds ; inds ; inds] ;
    st(inds3) = -[z(inds) ; z(inds) ; z(inds)] ;
    
    st = st-min(st(:)) ;
    st = st / max(st(:)) ;

    st = st.^(0.5) ;
    
    imagesc(st) ;
    saveas(h,sprintf('%s%d',filenames,i),'jpg') ;
end

end