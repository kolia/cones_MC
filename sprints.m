function s = sprints( format , varargin )

N = numel(varargin{1}) ;
for i=1:nargin-1
%     if iscell( varargin{i} ) ,  varargin{i} = cell2mat(varargin{i}) ; end
    varargin{i} = varargin{i}(:) ;
    if length(varargin{i}) ~= N
        error('sprints: numel(varargs{i}) must all be ==.') ;
    end
end

s = '' ;
format = ['%s' format] ;
arg = 's' ;
for i=1:nargin-1
    if iscell( varargin{i} )
        arg = sprintf('%s,varargin{%d}{j}',arg,i) ;
    else
        arg = sprintf('%s,varargin{%d}(j)',arg,i) ;
    end
end

for j=1:N
    s = eval(['sprintf(format,' arg ')']) ;
end

end