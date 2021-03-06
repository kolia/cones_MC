function gf = gaus_in_a_box(sigma)

% memo    = sparse([],[],[],1000,1000, ceil( (SS*3*sigma)^2 ) ) ;
gf      = @gib ;
 
    function out = gib(dx,dy)
    
%         out = zeros(size(in));  % preallocate output
%         [tf,loc] = ismember(in,x);  % find which in's already computed in x
%         ft = ~tf;  % ones to be computed
%         out(ft) = F(in(ft));  % get output values for ones not already in
%         % place new values in storage
%         x = [x in(ft(:).')];
%         y = [y reshape(out(ft),1,[])];
%         out(tf) = y(loc(tf));  % fill in the rest of the output values

        dx  = dx(:) ;
        dy  = dy(:) ;
        l   =  ones(size(dx)) ;
        O   = zeros(size(dx)) ;
        
        out = mvncdf([dx dy] + [l l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [O l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [l O],[],sigma.*[1 1]) ...
            + mvncdf([dx dy]        ,[],sigma.*[1 1]) ;
        out = out / sqrt(0.1083) ;
    end
end

% function gf = gaus_in_a_box(sigma,SS)
% 
% % memo    = sparse([],[],[],1000,1000, ceil( (SS*3*sigma)^2 ) ) ;
% memo    = cell(30,30,100,100) ;
% gf      = @gib ;
% 
%     function out = gib(dx,dy)
%         
% %         out = zeros(size(in));  % preallocate output
% %         [tf,loc] = ismember(in,x);  % find which in's already computed in x
% %         ft = ~tf;  % ones to be computed
% %         out(ft) = F(in(ft));  % get output values for ones not already in
% %         % place new values in storage
% %         x = [x in(ft(:).')];
% %         y = [y reshape(out(ft),1,[])];
% %         out(tf) = y(loc(tf));  % fill in the rest of the output values
% 
%     m = memo{100+dx*SS*2,100+dy*SS*2} ;
%     if isempty(m)
%         dx  = dx(:) ;
%         dy  = dy(:) ;
%         l   =  ones(size(dx)) ;
%         O   = zeros(size(dx)) ;
% 
%         out = mvncdf([dx dy] + [l l],[],sigma.*[1 1]) ...
%             - mvncdf([dx dy] + [O l],[],sigma.*[1 1]) ...
%             - mvncdf([dx dy] + [l O],[],sigma.*[1 1]) ...
%             + mvncdf([dx dy]        ,[],sigma.*[1 1]) ;
%         out = out / sqrt(0.1083) ;
%         memo{100+dx*SS*2,100+dy*SS*2} = out ;
%     else
%         out = m ;
%     end
%     end
% 
% end