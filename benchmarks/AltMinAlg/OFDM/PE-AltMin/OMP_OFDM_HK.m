% Modified version of OMP code for simulations.
% The algorithm is not changed. Only statistics collections are
% implemented. 
% Created by H.Kasai on July 01, 2018

function [ FRF, FBB, stats ] = OMP_OFDM_HK(Fopt, NRF, At)
    K = size(Fopt,3);
    FRF = [];
    Fres = Fopt;
    

    % set start time
    start_time = tic();
    
    for i = 1:NRF
        temp = 0;
        for k = 1:K
            PU(:,:,k) = At' * Fres(:,:,k);
            temp = temp + sum( abs(PU(:,:,k)).^2, 2 );
        end
        [aa,bb] = max(temp);
        FRF = [FRF , At(:,bb)];
        for k = 1:K
            FBB{k} = pinv(FRF) * Fopt(:,:,k);
            Fres(:,:,k) = (Fopt(:,:,k) - FRF * FBB{k}) / norm(Fopt(:,:,k) - FRF * FBB{k},'fro');
        end
    end
    
    % measure elapsed time
    elapsed_time = toc(start_time);  
    
    cost_val = 0;
    for k = 1:K      
        cost_val = cost_val + norm(Fopt(:,:,k) - FRF * FBB{k},'fro')^2;
    end 
    
    
    cost(1) = cost_val;
    time(1) = 0;      
    cost(2) = cost_val;
    time(2) = elapsed_time;
    
    stats.cost = cost;
    stats.time = time;        

end