function cone_map = greedy_cones( cone_map , speed )

warning off
if nargin<2
    speed = 0 ;
end

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
display_every = 1 ;
default( cone_map , 'plot_every'    , 0      )
default( cone_map , 'save_every'    , 500    )


%% PARAMETERS
% maxcones
  maxcones      = 3000 ;
  
fprintf('\n\nSTARTING greedy search')


% initializing variables
N_cones = 0 ;
X       = initialize_X(M0,M1,N_colors,SS,cone_map.cone_params.replusion_radii,...
                       cone_map.naive_LL,1,1) ;
X       = rmfield(X,'contact') ;
X.SS    = cone_map.cone_params.supersample ;
% X.LL_history = ones(ceil(M0*M1/10),1) ;

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% MAIN MCMC LOOP
t = cputime ;
tic
for jj=1:maxcones
    if speed
        [X,done] = greedy_hot(X,cone_map) ;
    else
        [X,done] = greedy(X,cone_map) ;
    end
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
    
    done = done | jj/(N_cones+1)>1.1 ;
    if ~mod(jj,save_every) || done
        to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
        try X = rmfield(X,'excluded') ; end
        X = remake_X(cone_map,X) ;
        try to_save.X = rmfield(X,{'invWW'}) ; end
%         try to_save.X = rmfield(X,{'WW'})    ; end
        save('result', 'to_save' )
        if done ,  break ;
        else clear to_save; end
    end
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.code           = file2str('greedy_cones.m') ;
cone_map.N_cones        = N_cones ;
% save('results','cone_map','X')

end