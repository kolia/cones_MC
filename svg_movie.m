function svg = svg_movie( init , dX , greedy, NICE )

data = 'data = [' ;
iterations = zeros(numel(dX),1) ;
for j=1:numel(iterations)
    iterations(j) = size(dX{j},1) ;
end

for i=[1 nnz(greedy):max(iterations)]
    idata = cell(numel(dX),4) ;
    for j=1:numel(dX)
        if size(dX{j},1)>=i
            for jj=1:3 , idata{j,jj} = [idata{j,jj} dX{j}(i,jj:3:end)] ; end
            idata{j,4} = [idata{j,4} repmat(j,[1 numel(idata{j,1})])] ;
            if  i == 1 , idata(j,:) = state2xyc( init{j} , j ) ; end
        end
    end
    if i==max(iterations), break ; end
    if i == 1
        ss = state2xyc(greedy,0) ;
        for k=1:4
            idata{k} = [idata{k} ss{k}] ; 
        end
    end
    ith = xycs2js( idata ) ;
%     if ~isempty( ith )
        data = [data '[' ith sprintf('],\n')] ;
%     end  
%     if ~mod(i,100), fprintf('%d\n',i) ; end
end
data = [data(1:end-1) ']'] ;

scale  = 200/max([size(NICE,1) size(NICE,2)]) ;
width  = min([200 200*size(NICE,2)/size(NICE,1)])+20 ;
height = min([200 200*size(NICE,1)/size(NICE,2)])+20 ;

evidence = print_evidence( NICE ) ;

fid = fopen('dance.js') ;
js = sprintf(fread(fid,'*char'),data) ;
fclose(fid) ;

fid = fopen('dancing_cones_source.svg') ;
svg = sprintf(fread(fid,'*char'),3*width,3*height,3*width,3*height,js,scale,scale,evidence) ;
fclose(fid) ;

% try
%     lsres = ls('dancing_cones.svg') ;
%     user_input = input('\nWarning: rm old dancing_cones.svg? (y/n): ','s') ;
%     switch user_input
%         case 'y'
%             delete('dancing_cones.svg')
%         otherwise
%             return
%     end
% end
fid = fopen('dancing_cones_movie.svg','w') ;
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