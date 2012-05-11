function [slowX,fastX] = swap_step(slowX,fastX,PROB,N_times,fastT)

% check_X(slowX)
% check_X(fastX)

% calculate overlaps of slowX cones on fastX exclusion disks
R = overlap_relation( fastX , slowX ) ;

% calculate at most N_times classes
classes  = equivalence_classes_direct(R,N_times) ;

if numel(classes)>0
    % propose swap moves N_times, choosing a different class every time
    for i=1:numel(classes)
        while 1
            ii = randi(numel(classes)) ;
            class = classes{ii} ;
            [proposed_slowX,proposed_fastX,class] = ...
                swap_closure( slowX, [1 1], fastX, fastT, PROB, class ) ;
            if ~isempty(proposed_slowX), break ; end
        end
        accept = metropolis_hastings(          slowX.ll+         fastX.ll,...
                                      proposed_slowX.ll+proposed_fastX.ll, 1) ;
%         if accept
%             fprintf('  swapdll: %.2f, slowdll: %.2f',...
%                 proposed_slowX.ll+proposed_fastX.ll-slowX.ll-fastX.ll,...
%                 proposed_slowX.ll-slowX.ll) ;  
%         end
        [slowX,fastX] = update_swap(slowX,fastX,proposed_slowX,proposed_fastX,accept) ;
        if accept, classes{ii} = class ; end
    end
end

end
