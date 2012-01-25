function cone_map = greedy_cones( cone_map )

% cone_map    = rmfield(cone_map,{'ROI'}) ;

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
plot_every      = 0 ;
display_every   = 1 ;

%% PARAMETERS
% maxcones
  maxcones      = 3000 ;
  
fprintf('\n\nSTARTING greedy search')


% initializing variables
N_cones = 0 ;
X       = initialize_X(M0,M1,N_colors,SS,cone_map.cone_params.replusion_radii,1,1) ;
X       = rmfield(X,'contact') ;
X.SS    = cone_map.cone_params.supersample ;

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% MAIN MCMC LOOP
t = cputime ;
tic
for jj=1:maxcones
    X = greedy(X,cone_map,@update_X) ;
    N_cones = numel(find(X.state>0)) ;
     
    % DISPLAY stdout
    if ~mod(jj,display_every)
        fprintf('\nCones:%4d, %4d  %.0f\tin %8.2f sec',jj,N_cones,X.ll,toc)
        tic
    end
    
    % DISPLAY plot
    if plot_every && ~mod(jj,plot_every)
        figure(h)
        plot_cones( X.state , cone_map ) ;
        title(sprintf('Iteration %d',jj),'FontSize',16)        
        drawnow
    end
    
    if N_cones<jj
        break ;
    end
%     save('results','cone_map')
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.code           = file2str('greedy_cones.m') ;
cone_map.N_cones        = N_cones ;
% save('results','cone_map','X')

end