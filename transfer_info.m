function B = transfer_info( A, B, byte_cutoff )
% Transfer all fields that are in a but not in b into b, only if their
% contents are smaller than byte_cutoff bytes (default is 10000).

if nargin<3,  byte_cutoff = 1e4 ; end

names = fieldnames(A) ;

for i=1:numel(names)
    name = names{i} ;
    if ~isfield(B,name)
        content = A.(name) ;
        content_whos = whos('content') ;
        content_size = content_whos.bytes ;
%         fprintf('%s: %d bytes\n',name,content_size)
        if content_size <= byte_cutoff
            B.(name) = content ;
        end
    end
end

end