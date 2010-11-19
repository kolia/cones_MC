function string = insert_string(s,filename,position)

fid = fopen(filename) ;

% Count backwards from end of file if position < 0
if position<0
    n_chars = numel( fread(fid, '*char') ) ;
    fclose(fid) ; fid = fopen(filename) ;
    position = n_chars + position ;
end

begin = fread(fid, position , '*char')' ;
ending = fread(fid, '*char')' ; fclose(fid) ;
string = [begin s ending] ;

end