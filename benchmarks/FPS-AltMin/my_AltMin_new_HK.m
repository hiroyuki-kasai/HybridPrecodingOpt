function [FRF, FBB, stats] = my_AltMin_new_HK( Fopt, C )
    [Nt, Ns] = size(Fopt);
    NRF = size(C,2);
    Nc = size(C,1)/NRF;
    [~,~,V] = svd(Fopt);
    FBB = [V';zeros(NRF-Ns,Ns)];

    mynorm = [Inf,0];
    
    index = 1;
    cost(index) = 0; % To Do : HK
    time(index) = 0;     
    
    % set start time
    start_time = tic();
    
    while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
        [alpha_k, value, s] = alpha_opt_new(real(Fopt*FBB'*C'));
        S = reshape(s,[Nt,Nc*NRF]);
        mynorm(1) = value + norm(imag(Fopt*FBB'*C'),'fro')^2;

        [U,~,V] = svds(alpha_k*Fopt'*S*C,NRF);
        FBB = V*U';
        mynorm(2) = norm(Fopt*FBB'*C' - alpha_k*S,'fro')^2;
        
        % measure elapsed time
        elapsed_time = toc(start_time);

        
        index = index + 1;
        cost(index) = real(mynorm(end));
        time(index) = elapsed_time;         
    end

    FBB = alpha_k*FBB;
    FRF = S*C;

    
    stats.cost = cost;
    stats.time = time;      
end