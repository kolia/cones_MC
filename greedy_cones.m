function cone_map = greedy_cones( GC_stas , cone_params , cone_map )
%% cone_map = greedy_cones( GC_stas , cone_params )
%  Do greedy search for cone positions.

addpath(genpath(pwd))

% size of region of interest
  [M0,M1,N_colors] = size(GC_stas(1).spatial) ;
  M2 = M0*M1 ;


%% PROBLEM SETUP: ganglion cell STAs

% supersample factor
SS  = cone_params.supersample ;

% set up Region of Interest
if isfield(cone_map,'ROI')
    ROI  = cone_map.ROI ;
    temp = cell(N_colors,1) ;
    for c=1:N_colors
        temp{c} = kron(ROI(:,:,c),ones(SS,SS)) ;
        temp{c} = temp{c}(:) ;
    end
    ROI  = find( cell2mat(temp) ) ;
else
    ROI  = 1:M0*M1*N_colors*SS^2 ;
end

% stereotyped cone receptive field
s = cone_params.sigma ;
r = cone_params.support_radius ;
cone_RF = exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2)' * exp(-0.5 * ((-SS*r:SS*r)/(SS*s)).^2) ;
cone_RF = cone_RF / sum(cone_RF(:)) ; clear r s

N    = M2*(SS^2)*N_colors ;

% Unpacking GC_stas into: STA, norms of STAs and N_spikes
N_GC = length(GC_stas) ;
STA_norm = zeros(N_GC,1) ;
N_spikes = zeros(N_GC,1) ;
STA      = zeros(N,N_GC) ;
for i=1:N_GC
    N_spikes(i) = length(GC_stas(i).spikes) ;
    temp = cell(N_colors,1) ;
    for c=1:N_colors
        temp{c} = kron(GC_stas(i).spatial(:,:,c),ones(SS,SS)) ;
        temp{c} = temp{c}(:) ;
    end
    STA(:,i)    = cell2mat(temp) ;
    STA_norm(i) = norm(STA(:,i)) ;
end

cell_consts = N_spikes ./ exp(STA_norm/2) * cone_params.stimulus_variance ;

%% SETUP for Log-LIKELIHOOD calculations

% file = sprintf('../resources/STA_W_super%d.mat',SS) ;

% try load(file) ; fprintf('\nloaded %s\n',file) 
% catch
    % W = matrix representing all M2 possible cone receptive fields
    STA_W = zeros(N,N_GC) ;
    fprintf('Calculating STA_W  --  counting up to %d:\n',M0*SS)
    for i=1:M0*SS
        fprintf('%d ',i)
        for j=1:M1*SS
            for c=1:N_colors
                BW      = zeros(M0*SS,M1*SS) ;
                BW(i,j) = 1 ;
                BW      = imfilter(BW,cone_RF) ;
                
                filter  = kron(cone_params.colors(c,:),BW(:)') ;
                
                STA_W((c-1)*M0*M1*SS^2 + (j-1)*M0*SS + i , :) = filter * STA ;
            end
        end
    end
    STA_W = STA_W .* repmat( (N_spikes ./ cell_consts)' ,N,1) ;
    
%     save(file,'STA_W')
% end

coneConv    = conv2(cone_RF,cone_RF) ;
colorDot    = cone_params.colors * cone_params.colors' ;
max_overlap = coneConv(ceil(size(coneConv,1)/2),ceil(size(coneConv,2)/2)) ;

prior_LL = @(X)-1e12*sum(sum( triu(X.overlaps,1) > max_overlap*0.1 )) ;


%% GREEDY SOLUTION
GREED.state  = zeros( 1 , N ) ;
best_LL      = -Inf ;

flip_LL      = @(X,flips)flip_color_LL( ...
    X , flips , prior_LL , cell_consts , STA_W' , coneConv , colorDot , [M0 M1]*SS) ;


fprintf('\nGREEDY cone finding:\n')
tic
while 1
    GREED = greedy( GREED , flip_LL , ROI ) ;
    fprintf('\nCONES:%2d \t \t increase in LL:%f',sum(GREED.state),GREED.ll-best_LL)
    if GREED.ll<=best_LL
        break
    else
        best_LL = GREED.ll ;
    end
%     find(GREED.state)
%     GREED.overlaps - diag(diag(GREED.overlaps))
end
fprintf('\nGREEDY SOLUTION found in %.1f sec\n',toc)

cone_map.stas        = GC_stas ;
cone_map.cone_RF     = cone_RF ;
cone_map.cone_params = cone_params ;
cone_map.prior_LL    = func2str(prior_LL) ;
cone_map.max_overlap = max_overlap ;
cone_map.GREED       = GREED ;
cone_map.coneConv    = coneConv ;
cone_map.colorDot    = colorDot ;
cone_map.ROI_super   = ROI ;
end