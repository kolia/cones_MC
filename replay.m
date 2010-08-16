function replay(nums,skip,speed)
% speed in frames/sec

if nargin<1 ,   nums  = [1;2;3;4;5;6;7;8]   ; end
if nargin<2 ,   skip  = 10                  ; end
if nargin<3 ,   speed = 10                  ; end

load('QC.mat')
load('gcm.mat')
[gx,gy] = find(gcm.X.state) ;

M0 = 26 ;
M1 = 46 ;
SS = 4  ;

NN  = numel(nums) ;
dX  = cell(NN,1) ;

n = 0 ;
N = 0 ;
states = cell(NN,1) ;
for i=1:NN
    stats = load(sprintf('stats_%d.mat',nums(i))) ;
    dX{i} = stats.results{1}.dX ;
    clear stats
    n = max(n,numel(max(dX{i}(:,1:3:end)))) ;
    N = max(N,find(max(dX{i},[],2),1,'last')) ;
    states{i} = zeros(M0*SS,M1*SS) ;
end

scrsz = get(0,'ScreenSize');
hf = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;

colormap('pink')
imagesc(QC) ;
hold on

plot(gy,gx,'w+','MarkerSize',7)  

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
    
    if i>4500  && ~mod(i,skip)
        tic ;

        for ii=1:NN
            s = symbols{ii} ;
            for cc=1:3
                c = colors{cc} ;
                
                [ix,iy] = find(states{ii} == cc) ;
                h{cc,ii} = plot(iy,ix,sprintf('%s%s',c,s),'MarkerSize',5) ;
            end
        end
        
        title(sprintf('Iteration %d',i),'FontSize',16)
        drawnow expose

        pause(1/speed - toc) ;
        
        if i<N-skip ,  delete(cell2mat(h(:))) ; end        
        
    end    
end