function [FRF, FBB, stats] = OMP_NB_HK( Fopt, NRF, At )

    FRF = [];
    Fres = Fopt;
    
    % set start time
    start_time = tic();
    
    for k = 1:NRF
        PU = At' * Fres;
    %     [aa,bb] = max(diag(PU * PU'));
        [aa,bb] = max(sum( abs(PU).^2, 2 ));
        FRF = [FRF , At(:,bb)];
        FBB = pinv(FRF) * Fopt; %use pseudoinverse to avoid the inverse of a possible singular matrix
        Fres = (Fopt - FRF * FBB) / norm(Fopt - FRF * FBB,'fro');
    end
    
    % measure elapsed time
    elapsed_time = toc(start_time);  
    
    cost_val = norm(Fopt(:,:) - FRF * FBB,'fro')^2;    
    
    cost(1) = cost_val;
    time(1) = 0;      
    cost(2) = cost_val;
    time(2) = elapsed_time;
    
    stats.cost = cost;
    stats.time = time;       

end

