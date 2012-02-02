function stas = restrict_ROI( stas, xrange, yrange )
% Restrict stas to new Region Of Interest [x X], [y Y].
% Only keep those cells for which the sta extremum is within the new ROI.

keep = zeros(numel(stas),1) ;
for c=1:numel(stas)
    [~,I] = max( sum(abs(stas(c).spatial),3) ) ;
    if xrange(1)<=I(1) && I(1)<=xrange(2) && yrange(1)<=I(2) && I(2)<=yrange(2)
        keep(c) = 1 ;
    end
end
stas = stas(logical(keep)) ;

for c=1:numel(stas)
    stas(c).spatial = stas(c).spatial(xrange(1):xrange(2),yrange(1):yrange(2),:) ;
end

end