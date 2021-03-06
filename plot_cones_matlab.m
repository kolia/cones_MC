function h = plot_cones_matlab( states , cone_map )
% matlab plot of a cone configuration, using states and cone_map.NICE

if isnumeric(states)
    states = { states } ;
end

if iscell(states) && isfield(states{1},'state')
    for i=1:numel(states)
        states{i} = states{i}.state ;
    end
end

% if issparse(states{1})
%     for ii=1:numel(states)
%         states{ii} = states{ii}( states{ii}>0 ) ;
%     end
% end

symbols = {'o' ; 's' ; '*' ; 'x' ; 'p' ; 'h' ; 'd' ; '.'} ;
colors  = {'r' ; 'g' ; 'b' } ;

NN = numel(states) ;
h  = cell(3,NN) ;

close all
figure('visible','off')

if nargin>1 && isfield( cone_map ,'NICE' )
    colormap('pink')
    imagesc( cone_map.NICE.^(0.6) ) ;
end

roughsize = max(size(cone_map.NICE)) ;

hold on
for ii=1:numel(states)
    s = symbols{ii} ;
    for cc=1:3
        c = colors{cc} ;
        
        [ix,iy] = find(states{ii} == cc) ;
        h{cc,ii} = plot(iy,ix,sprintf('%s%s',c,s),'MarkerSize',ceil(6*640/roughsize)) ;
    end
end

if nargin>1
    if isfield(cone_map,'greedy')
        cone_map = cone_map.greedy ;
    end
    if isnumeric(cone_map)
        [gx,gy] = find( cone_map ) ;
        plot(gy,gx,'w+','MarkerSize',7) ;
    end
end

end