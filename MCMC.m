function bestX = MCMC( cone_map , ID )

if nargin>1 ,   cone_map.ID = ID ;     end

M0          = cone_map.M0 ;
M1          = cone_map.M1 ;
cone_map.SS = cone_map.cone_params.supersample ;
SS          = cone_map.SS ;
N_colors    = cone_map.N_colors ;

initX = initialize_X( M0, M1, N_colors, SS, ...
                        cone_map.cone_params.replusion_radii, 1, 1) ;

default( cone_map , 'N_iterations' , 100000) ;
default( cone_map , 'X' , initX   )
default( cone_map , 'q' , 0.999 )

default( cone_map , 'plot_every'    , 0     )
default( cone_map , 'plot_skip'     , 100   )
default( cone_map , 'display_every' , 100   )
default( cone_map , 'save_every'    , 200   )
default( cone_map , 'ID'            , 0     )
default( cone_map , 'N_best'        , 1     )
default( cone_map , 'max_time'      , 1000  )


% sparse int matrix, with number of out-of-border adjacencies
cone_map.outofbounds = sparse([],[],[],M0*SS,M1*SS,2*(M0+M1)*SS) ;
cone_map.outofbounds(:,[1 M1*SS]) = 1 ;
cone_map.outofbounds([1 M0*SS],:) = cone_map.outofbounds([1 M0*SS],:) + 1 ;


% initializing variables
bestX = cell(N_best,1) ;
for i=1:N_best
    bestX{i} = X ;
end

if plot_every
scrsz = get(0,'ScreenSize');
h = figure('Position',[1 scrsz(4)*0.7*0.5 1500*0.5 1200*0.5]) ;
end

% MAIN MCMC LOOP
fprintf('\n\nStochastic Tunneling with alpha = %g  :\n',aleph)
t = cputime ;
tic

jj = 0 ;
while 1

    if mean_f<0
        if - mean_f < 0.1
            beta = beta/0.96 + 0.01 ;
        else
            beta = beta*0.9 ;
        end
    end
    
    % reinitialize if stuck
    if beta>1e20
        X       = initX ;
        beta    = init_beta ;
    end
    
    [ d, X ] = flip_MCMC( struct, X, move( X , 2 , q , cone_map ),...
                          @update_X , @(trial)get_ll(trial,beta) ) ;
    n_cones = numel(find(X.state>0)) ;

    N_mean = 50 ;
    mean_f = ( mean_f * (N_mean-1) - (X.ll - bestX{1}.ll) )/N_mean ;    
    
    if N_best
        for b=1:N_best
            if X.ll>bestX{b}.ll && ...
                    ( b==1 || X.ll<bestX{b-1}.ll)
                for bb=b+1:N_best
                    bestX{bb} = bestX{bb-1} ;
                end
                bestX{b} = X ;
            end
        end
    end
    
    % DISPLAY stdout
    if ~mod(jj,display_every)
        fprintf('Iteration:%4d of %d  %4d cones    %6.0f L   %6.0f best   %10g beta   %10g dLL   %8.2f sec\n',...
                            jj,N_iterations,n_cones,X.ll,...
                            bestX{1}.ll,beta,mean_f,toc)
        tic
    end
    
    % DISPLAY plot
    if ~mod(jj,plot_every)
        figure(h)
        plot_cones( X.state , cone_map ) ;
        title( sprintf('After %d Stochastic Tunneling iterations',jj),...
               'FontSize' , 24 )
        % set(get(gca,'Title'),'Visible','on')
        drawnow
    end
    
    if ~mod(jj,save_every)
        save(sprintf('bestX_%d',ID), 'bestX')
    end
    
    jj = jj + 1 ;

    if jj>N_iterations || cputime-t>max_time ,  break ;  end
end
fprintf('\ndone in %.1f sec\n\n', cputime - t) ;

cone_map.X              = X ;
cone_map.bestX          = bestX ;
cone_map.code.tunnel    = file2str('tunnel.m') ;
save(sprintf('bestX_%d',ID), 'bestX')

end


function default( s , name , value )

if isfield( s , name )
    assignin( 'caller' , name , s.(name) ) ;
else
    assignin( 'caller' , name , value ) ;
end

end