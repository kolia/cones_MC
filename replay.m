function replay(num,speed,skip)
% speed in iterations/sec

if nargin<1 ,   num = 1 ;     end
if nargin<2 ,   speed = 100 ; end
if nargin<3 ,   skip  = 10  ; end

load('/Users/kolia/Desktop/QC.mat')
load(sprintf('/Users/kolia/Desktop/stats_%d.mat',num))
load('/Users/kolia/Documents/github/cones_MC/gcm.mat')

M0 = 26 ;
M1 = 46 ;
SS = 4  ;

dX = results{1}.dX ;
clear results swap_stas

% GQC = QC ;
% for c=1:3
%     [ix,iy] = find(gcm01.X.state == c) ;
%     for cc=1:3
%         GQC(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(cc-1)) = 0 ;
%     end
%     GQC(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(c-1)) = 1 ;
% end

[gx,gy] = find(gcm.X.state) ;

scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;

figure(h)

state = zeros(M0*SS,M1*SS) ;

n     = numel(max(dX(:,1:3:end))) ;

for ii=1:size(dX,1)

    i = 1+mod(ii-1,7000) ;
    if i == 1
        state = zeros(M0*SS,M1*SS) ;
    end
    
    tic ;
    for j=0:n-1
        if dX(i,1+3*j)
            state(dX(i,1+3*j),dX(i,2+3*j)) = dX(i,3+3*j) ;
        end
    end
    
    if ~mod(i,skip)
        colormap('pink')
        GGG = QC ;
        for c=1:3
            [ix,iy] = find(state == c) ;
            for cc=1:3
                GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(cc-1)) = 0 ;
            end
            GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(c-1)) = 1 ;
        end
        
        imagesc(GGG) ;
        
        hold on
        plot(gy,gx,'w.','MarkerSize',5)
        
        title(sprintf('Iteration %d',i),'FontSize',16)
        drawnow

        if ii>=3500
            waitforbuttonpress
        end
    
    end    
    pause( max( 1/speed - toc , 0 ) )
end