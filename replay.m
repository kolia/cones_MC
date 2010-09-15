function replay( dX, greedy, cone_map   , skip, speed )

if nargin<4 ,   skip  = 10  ; end  % iterations/frame
if nargin<5 ,   speed = 10  ; end  % frames/sec

[gx,gy] = find(greedy) ;

M0 = cone_map.M0 * cone_map.cone_params.supersample ;
M1 = cone_map.M1 * cone_map.cone_params.supersample ;
NN = numel(dX) ;

n = 0 ;
N = 0 ;
states = cell(NN,1) ;
for i=1:NN
    n = max(n,numel(max(dX{i}(:,1:3:end)))) ;
    N = max(N,find(max(dX{i},[],2),1,'last')) ;
    states{i} = zeros(M0,M1) ;
end

scrsz = get(0,'ScreenSize');
hf = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;

colormap('pink')
imagesc( cone_map.NICE ) ;
hold on

symbols = {'.' ; 's' ; '*' ; 'x' ; 'p' ; 'h' ; 'd' ; 'o'} ;
colors  = {'r' ; 'g' ; 'b' } ;

h = cell(3,NN) ;

for i=1:N
    for j=0:n-1
        for ii=1:NN
            if dX{ii}(i,1+3*j)
                states{ii}(dX{ii}(i,1+3*j),dX{ii}(i,2+3*j)) = dX{ii}(i,3+3*j) ;
            end
        end
    end
    
    if ~mod(i,skip)
        tic ;

        figure(hf)
        
        for ii=1:NN
            s = symbols{ii} ;
            for cc=1:3
                c = colors{cc} ;
                
                [ix,iy] = find(states{ii} == cc) ;
                h{cc,ii} = plot(iy,ix,sprintf('%s%s',c,s),'MarkerSize',5) ;
            end
        end        
        plot(gy,gx,'w+','MarkerSize',7)
        
        title(sprintf('Iteration %d',i),'FontSize',16)
        drawnow expose
        pause(1/speed - toc) ;        
        if i<N-skip ,  delete(cell2mat(h(:))) ; end        
        
    end    
end