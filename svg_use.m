function svg = svg_use( x , y , id )

if iscell(x) ,  x = cell2mat(x) ; end
if iscell(y) ,  x = cell2mat(y) ; end
    
x = x(:) ;
y = y(:) ;
N = length(x) ;
if length(y) ~= N, error('svg_use: numel(args) must be ==.') ; end

svg = '' ;
for i=1:min([numel(x) numel(y)])
    svg = sprintf('%s<use xlink:href="#%d" transform="translate(%f %f)"/>\n',...
                   svg,                id(i),                 x(i), y(i) ) ;
end

end