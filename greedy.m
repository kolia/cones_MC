function cone_map = GREEDY( cone_map , speed )

warning off
if nargin<2     % default is regular greedy (anything but 0: 'hot' version)
    speed = 0 ;
end

% plot, display, and save every N iterations (0 = never)
display_every = 1 ;
default( cone_map , 'plot_every'    , 0      )
default( cone_map , 'save_every'    , 500    )

% no need to track_contacts, greedy does not shift cones, it just adds them
if ~isfield(cone_map , 'track_contacts'), cone_map.track_contacts = false ; end

if speed
    fprintf('\n\nSTARTING ''hot'' greedy search')
else
    fprintf('\n\nSTARTING greedy search')
end

% initializing variables
X = cone_map.initX ;

% if not tracking contacts, contact field is not needed
if ~cone_map.track_contacts, try X   = rmfield(X,'contact') ; end ; end

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% GREEDY addition of cones, one at a time
t = cputime ;
tic
for jj=1:cone_map.initX.maxcones

    % try adding a cone
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
    
    % SAVE to disk
    save_now = done | jj/(N_cones+1)>1.1 ;
    if ~mod(jj,save_every) || save_now
        % use less disk space: remove large data structures
        to_save = rmfield(cone_map,{'STA','initX','sparse_struct'}) ;
        try X = rmfield(X,'excluded') ; end
        X = remake_X(cone_map,X) ;
        try to_save.X = rmfield(X,{'invWW'}) ; end
        save('result', 'to_save' )
        if save_now ,  break ;
        else clear to_save; end
    end
end
fprintf('\ndone in %.1f sec\n\n',cputime - t) ;

cone_map.X              = X ;
cone_map.code           = file2str('greedy_cones.m') ;
cone_map.N_cones        = N_cones ;
% save('results','cone_map','X')

end