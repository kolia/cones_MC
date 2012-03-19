function save_svg_plot(svg,filename,batik)

if nargin<3
    batik = '../batik/' ;
end

fid = fopen(filename,'w') ;
fwrite(fid,svg) ; fclose(fid) ;
system(['java -jar ' batik 'batik-rasterizer.jar -m application/pdf ' filename])

end