function make_sta_plots( sta , invww, filenames )
% Given a cell struct of stas, for example as output by denoised_sta, save
% imagesc plots to current directory.  Negative regions in STA are shown in
% white.

NGC = size( sta, 2 ) ;

h = figure ;

for i=1:NGC
    st = plotable_evidence( sta{1,i} ) ;    
    st = st.^(0.5) ;

    try
        imagesc(st) ;
        filename = sprintf('%s_%d_STA',filenames,i) ;
        saveas(h,filename,'jpg') ;
    catch E
        fprintf('Error while attempting to plot %s',filename)
        disp(E)
    end
    
    try
        std = sqrt( sta{2,i} - sta{1,i}.^2 + invww ) ;
        std = std-min(std(:)) ;
        std = std / max(std(:)) ;
        imagesc(std)
        filename = sprintf('%s_%d_std',filenames,i) ;
        saveas(h,filename,'jpg') ;
    catch E
        fprintf('Error while attempting to plot %s',filename)
        disp(E)
    end
end