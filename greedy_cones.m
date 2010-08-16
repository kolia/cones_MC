function cone_map = greedy_cones( cone_map )

LL          = cone_map.LL ;
cone_map    = rmfield(cone_map,{'LL' 'ROI'}) ;

M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
cone_map.SS = cone_map.cone_params.supersample ;
SS          = cone_map.SS ;
N_colors    = cone_map.N_colors ;

% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;

%% plot and display info every ?
plot_every      = 1 ;
display_every   = 1 ;

%% PARAMETERS
% D             exclusion distance
  D             = 9.2 ;
% maxcones
  maxcones      = 150 ;

% FACTOR        arbitrary factor for W
  FACTOR        = 0.1 ;
  
cone_map.STA_W      = cone_map.STA_W    * FACTOR    ;
cone_map.coneConv   = cone_map.coneConv * FACTOR^2  ;

  
fprintf('\n\nSTARTING greedy search')

% initializing variables
n_cones = 0 ;
X       = initialize_X(M0,M1,N_colors,SS,D,1,maxcones) ;
X       = rmfield(X,'contact') ;
X.SS    = cone_map.cone_params.supersample ;
    
results.iteration = 0 ;
results.dX        = sparse([],[],[],maxcones,3*maxcones) ;


if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

GG0 = LL - min(LL(:)) ;
GG0 = ( GG0 / max(GG0(:))) .^ 0.7 ; % !!! enhance low values visually


% MAIN MCMC LOOP
t = cputime ;
tic
for jj=1:maxcones
    [results,X] = greedy(results,X,cone_map,@update_X) ;
    n_cones = numel(find(X.state>0)) ;
     
    % DISPLAY stdout
    if ~mod(jj,display_every)
        fprintf('\nCones:%4d , %4d  \tin %8.2f sec',jj,n_cones,toc)
        tic
    end
    
    % DISPLAY plot
    if ~mod(jj,plot_every)
        figure(h)
        colormap('pink')
        GGG = GG0 ;
        for c=1:3
            [ix,iy] = find(X.state == c) ;
            for cc=1:3
                GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(cc-1)) = 0 ;
            end
            GGG(ix + M0*SS*(iy-1) + M0*SS*M1*SS*(c-1)) = 1 ;
        end
        imagesc(GGG) ;
        title(sprintf('Iteration %d',jj),'FontSize',16)
        drawnow
    end
    
    if n_cones<jj
        break ;
    end
    save('results_f01','results','cone_map')
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.code           = file2str('greedy_cones.m') ;
cone_map.n_cones        = n_cones ;
save('results','results','cone_map')

end