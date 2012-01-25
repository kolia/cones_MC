function gf = gaus_in_a_box_memo(sigma,SS,supprt)

% memo    = sparse([],[],[],1000,1000, ceil( (SS*3*sigma)^2 ) ) ;
memo    = cell(ceil(supprt+2),ceil(supprt+2)) ;
gf      = @gib ;

    function out = gib(i,j)
    di = rem(i,1) ;
    dj = rem(j,1) ;

    out = memo{0.5+di*SS,0.5+dj*SS} ;
    if isempty(out)
        
        ox    = 1+(floor(i-supprt):floor(i+supprt)) ;
        oy    = 1+(floor(j-supprt):floor(j+supprt)) ;
        dx    = repmat(ox(:),numel(oy),1) ;
        dy    = reshape( repmat(oy,numel(ox),1) , [] , 1 ) ;

        dx  = i-dx(:) ;
        dy  = j-dy(:) ;
        l   =  ones(size(dx)) ;
        O   = zeros(size(dx)) ;

        out = mvncdf([dx dy] + [l l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [O l],[],sigma.*[1 1]) ...
            - mvncdf([dx dy] + [l O],[],sigma.*[1 1]) ...
            + mvncdf([dx dy]        ,[],sigma.*[1 1]) ;
        out = out / sqrt(0.1083) ;
        memo{0.5+di*SS,0.5+dj*SS} = out ;
    end 
    end
end