load george/full_stas
load george/time_courses
load george/datarun

rf = zeros( 160, 160, 3 ) ; %0 * squeeze(full_stas{1}(:,:,:,1)) ;

k = 0 ;
r = struct ;
for i=1:376
    if numel(full_stas{i})>0
        k = k+1 ;
        r(k).spatial = rf ;
        r(k).spikes  = datarun.spikes{i} ;
        r(k).cell_index = i ;
        for c=1:3
            timecourse = time_courses{i}(:,c) ./ sqrt( sum(time_courses{i}(:,c).^2) ) ;
            for j=1:7
                r(k).spatial(:,:,c) = r(k).spatial(:,:,c) + ...
                         squeeze(full_stas{i}(81:240,61:220,c,j))*timecourse(j) ;
            end
        end
    end
end

stas = r ;
save('george/stas','stas')