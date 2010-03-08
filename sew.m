function quilts = sew( r , accessors )

if nargin<2
    if isfield(r{1},'GREED')
        accessors = {@(rr)rr.GREED.state} ;
    elseif isfield(r{1},'accumulated')
        accessors = {@(rr)rr.accumulated{1}(3:end)} ;
    end
end

N        = length(r) ;
M        = length(accessors) ;
S        = r{1}.cone_params.supersample ;
quilts   = cell(M,1) ;
position = cell(N,1) ;

sizes = r{1}.cut.total_sizes ;
sizes(1:2) = sizes(1:2) * r{1}.cone_params.supersample ;

% overX = ceil( (r{1}.position(1,end) - r{2}.position(1,1) + 1) * S/2 ) ;
% overY = 12 ;

sx    = ( r{1}.cut.position{1}(end) - r{1}.cut.position{1}(1) + 1 ) ;
sy    = ( r{1}.cut.position{2}(end) - r{1}.cut.position{2}(1) + 1 ) ;

overX = floor((sx - r{1}.cut.X(2) + 1)*S/2) ;
overY = floor((sy - r{1}.cut.Y(2) + 1)*S/2) ;

split = [zeros(N,1) ones(N,1)*sizes(1)+1 zeros(N,1) ones(N,1)*sizes(2)+1] ;
for i=1:N
    position{i} = cell(2,1) ;
    for k=1:2
        position{i}{k} = r{i}.cut.position{k}(:) ;
        position{i}{k} = (position{i}{k}(1)-1)*S+1:position{i}{k}(end)*S ;
    end
end

nX = position{1}{1}(end) - position{1}{1}(1) + 1 ;
nY = position{1}{2}(end) - position{1}{2}(1) + 1 ;

X = (r{1}.cut.X-1) * S + 1 ;
Y = (r{1}.cut.Y-1) * S + 1 ;

for j=1:length(Y)
    for i=1:length(X)
        ii = (j-1)*length(X)+i ;
        if i > 1
            split(ii,1) = floor( (X(i-1) + nX - 1 + X(i))/2 ) ;
        end
        if i < length(X)
            split(ii,2) = floor( (X(i) + nX - 1 + X(i+1))/2 ) ;
        end
        if j > 1
            split(ii,3) = floor( (Y(j-1) + nY - 1 + Y(j))/2 ) ;
        end
        if j < length(Y)
            split(ii,4) = floor( (Y(j) + nY - 1 + Y(j+1))/2 ) ;
        end
    end
end


for j=1:M
    quilts{j} = zeros(sizes) ;
%     counts    = zeros(sizes) ;
        
    for i=1:N
%         quilts{j} = quilts{j} + reshape( accessors{j}(r{i}) , sizes ) ;
        
%         n = length(position{1}) ;
%         m = length(position{2}) ;
        
%         size( accessors{j}(r{i}) )
%         [n m 3]
        
%         A  = reshape( accessors{j}(r{i}) , [n m 3] ) ;
        A  = reshape( accessors{j}(r{i}) , sizes ) ;

%         subplot(2,1,1)
%         imagesc( A ) 

        A(1:split(i,1)  ,:,:) = 0 ;
        A(split(i,2)+1:end,:,:) = 0 ;
        A(:,1:split(i,3)  ,:) = 0 ;
        A(:,split(i,4)+1:end,:) = 0 ;
        A = A(1:sizes(1),1:sizes(2),:) ;
        
%         [position{i}{1} position{i}{2} split(i,:)]
        
%         waitforbuttonpress
%         close all
        
        quilts{j} = quilts{j} + A ;

%         subplot(2,1,2)
%         imagesc(quilts{j} / max(quilts{j}(:)))
%         waitforbuttonpress
        
%         quilts{j}(position{1},position{2},1:3) = ...
%             quilts{j}(position{1},position{2},1:3) + A ;
        
%         counts(position{1},position{2},1:3) = ...
%             counts(position{1},position{2},1:3) + 1 ;
    end
    
%     quilts{j} = quilts{j}>0 ;
end