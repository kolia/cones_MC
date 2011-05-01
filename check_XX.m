function check_XX(X)

if numel(find(X.state>0)) ~= size(X.invWW,1)
   'oups'
end

end