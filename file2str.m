function s = file2str(filename)

fid = fopen(filename,'r') ;
s   = fread(fid,'*char')' ;
fclose(fid) ;

end