function masks = make_masks(D)

s = floor(2*D + 5) ;
c = floor(D   + 3) ;

disks = cell(5,1) ;
for d=1:5
    disks{d} = zeros(s) ;
end

shift   = cell(5,1) ;
shift{1}= [ 0  0] ;
shift{2}= [ 1  0] ;
shift{3}= [ 0  1] ;
shift{4}= -shift{2} ;
shift{5}= -shift{3} ;

for i=1:s
    for j=1:s
        for d=1:5
            if sqrt( (i-c+shift{d}(1))^2 + (j-c+shift{d}(2))^2 ) <= D
                disks{d}(i,j) = 1 ;
            end
        end
    end
end

for d=1:4
    [x,y] = find( disks{1+d} & ~disks{1} ) ;
    masks.cardinal{d} = [x y] - c ;
end

% % testing with plot
% im = zeros(s,s,3) ;
% for d=1:3
%     inds = s*s*(d-1)  + masks.cardinal{d}(:,1)+c + s * (masks.cardinal{d}(:,2)+c-1) ;    
%     im( inds ) = 1 ;
% end
% imagesc(im)