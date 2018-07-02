% Modified version of PE_AltMin code for simulations.
% The algorithm is not changed. Only statistics collections are
% implemented. 
% Created by H.Kasai on July 01, 2018

function [ FRF,FBB, stats ] = PE_AltMin_HK( Fopt, NRF, FRF_in )

    [Nt, Ns, K] = size(Fopt);
    mynorm = [];
    %FRF = exp( 1i * unifrnd (0,2*pi,Nt,NRF) );
    % FBB = zeros(NRF, Ns, K);

    %FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    if FRF_in == 0
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    else
        FRF = FRF_in;
    end 

    init_cost = 0;
    for k = 1:K      
        FBB(:,:,k) = pinv(FRF) * Fopt(:,:,k); 
        init_cost = init_cost + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
    end
     
    index = 1;
    cost(index) = init_cost;
    time(index) = 0;    

    % set start time
    start_time = tic();
    
    while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-1)
        mynorm = [0,0];
        temp = zeros(Nt, NRF);
        for k = 1:K
            [U,S,V] = svd(Fopt(:,:,k)'*FRF);
            FBB(:,:,k) = V(:,[1:Ns])*U';
            mynorm(1) = mynorm(1) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
            %mynorm(1) = mynorm(1) + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
            temp = temp + Fopt(:,:,k) * FBB(:,:,k)';
        end

        FRF = exp(1i * angle(temp));
        
        % measure elapsed time
        elapsed_time = toc(start_time);
        
        for k = 1:K
            mynorm(2) = mynorm(2) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
            %mynorm(2) = mynorm(2) + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
        end
        
        index = index + 1;
        cost(index) = real(mynorm(2));
        time(index) = elapsed_time;        
    end
    
    
    stats.cost = cost;
    stats.time = time;    
end