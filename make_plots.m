function make_plots( result, varargin )
%    make_plots( X )
% or make_plots( X, cone_map )
% or make_plots( result )  
%        where X = result.X and cone_map = result
%        as returned by GREEDY, MCMC, and CAST
%
% makes whatever plots it can, saves them in folder
% ['plots_for_' cone_map_string(result)]

if iscell(result)
    result = result{1} ;
end

bestX = [] ;
if isfield(result,'bestX')
    bestX = result.bestX ;
end

if isfield(result,'state') % result is an X
    X = result ;
    if ~isempty(varargin)
        cone_map = result ;        
    else
        cone_map = remake_cone_map( X ) ;
    end
else                       % result is a cone_map
    X = result.X ;
    if iscell(X), X = X{1} ; end
    cone_map = rmfield(result,'X') ;
end
clear result

folder_name = cone_map_string(cone_map);

% place all plots in a folder
filename = ['plots_for_' folder_name] ;
mkdir( filename )
here = pwd ;
cd(filename)

try
    % if possible, make cone field plot and movie; this requires dX
    % CAST runs have a dX field that records all moves
    if isfield( X, 'dX' )
        try load('./confident')
        catch
            cone_map = remake_initX( cone_map, true ) ;
            fprintf('making confident.mat\n')
            % select which iterations to use; 10000 burnin, then every 20 iters
            selector = @(n) (n>1000) && (mod(n,20) == 0) ;
            confident = confident_cones( cone_map.initX , X.dX , cone_map , selector ) ;
        end
        save('confident','confident') ;
            
        svg_movie( {cone_map.initX.state} , {X.dX} , cone_map.initX.state, cone_map.NICE ) ;
        
        plot_cone_field( confident , cone_map )
    end

    if ~isempty(bestX)
        for k=1:numel(bestX)
            plot_cones_matlab( bestX{k}.state, cone_map ) ;
            saveas(gcf,sprintf('best_cone_configuration_matlab_%d',k),'jpg')
            plot_cone_field( bestX{k} , cone_map, sprintf('best_cone_configuration_%d',k) )
        end
    else
        plot_cones_matlab( X.state, cone_map ) ;
        saveas(gcf,'final_cone_configuration_matlab','jpg')
        plot_cone_field( X , cone_map, 'final_cone_configuration' )
    end

    % % [sta,invww] = denoised_sta( X.initX , X.dX , cone_map, selector ) ;
    % % make_sta_plots( sta , invww, 'denoised' )

    cd(here)
catch err
    cd(here)
    rethrow(err);
end

end