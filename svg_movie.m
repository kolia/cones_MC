function svg = svg_movie( init , dX , greedy )

data = 'data = [' ;
for i=1:1e10
    idata = cell(numel(dX),4) ; run = 0 ;
    for j=1:numel(dX)
        if size(dX{j},1)>=i
            if  i == 1 , idata(j,:) = state2xyc( init{j} , j ) ; end
            for jj=1:3 , idata{j,jj} = [idata{j,jj} dX{j}(i,jj:3:end)] ; end
            idata{j,4} = [idata{j,4} repmat(j,[1 numel(idata{j,1})])] ;
            run = 1 ;
        end
    end
    if ~run ,  break ; end
    if i == 1 , idata = [idata ; state2xyc(greedy,0)] ; end
    
    ith = xycs2js( idata ) ;
    if ~isempty( ith )
        data = [data '[' ith sprintf('],\n')] ;
    end
end
data = [data(1:end-1) ']'] ;

svg = insert_string( data , which('dance.js') , 356 ) ;
svg = insert_string( svg  , which('dancing_cones_source.svg') , 172 ) ;

try
    lsres = ls('dancing_cones.svg') ;
    user_input = input('\nWarning: rm old dancing_cones.svg? (y/n):','s') ;
    switch user_input
        case 'y'
            delete('dancing_cones.svg')
        otherwise
            return
    end
end
fid = fopen('dancing_cones.svg','w') ;
fwrite(fid,svg) ; fclose(fid) ;

end


function xyc = state2xyc(state,ii)
[x,y,c] = find( state ) ; xyc = {x(:)' y(:)' c(:)' repmat(ii,[1 numel(x)])} ;
end


function js = xycs2js( xycs )
js = '' ;
for i=1:size(xycs,1)
    inds = find(xycs{i,1}) ;
    for j=inds
        js = [js '[ ' num2str(xycs{i,2}(j)) ' , ' num2str(xycs{i,1}(j)) ' , ' ...
                      num2str(xycs{i,3}(j)) ' , ' num2str(xycs{i,4}(j)) ']  ,  '] ;
    end
end
if ~isempty(js) ,  js = js(1:end-5) ; end
end