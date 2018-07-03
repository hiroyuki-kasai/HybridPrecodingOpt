% Modified version of PE_AltMin code for simulations.
% The algorithm is not changed. Only statistics collections are
% implemented. 
% Created by H.Kasai on July 01, 2018

function [ FRF,FBB, stats ] = PE_AltMin_NB_HK( Fopt, NRF, FRF_in )

    [Nt, Ns] = size(Fopt);

    if FRF_in == 0
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    else
        FRF = FRF_in;
    end 

    FBB(:,:) = pinv(FRF) * Fopt(:,:); 
    init_cost = norm(Fopt(:,:) - FRF * FBB(:,:),'fro')^2;

    index = 1;
    cost(index) = init_cost;
    time(index) = 0; 
    
    mynorm = [];

    % set start time
    start_time = tic();
    
    while (isempty(mynorm) || abs( mynorm(end) - mynorm(end-1) ) > 1e-3)
        
        temp = zeros(Nt, NRF);
        [U,S,V] = svd(Fopt'*FRF);
        FBB = V(:,[1:Ns])*U';
        mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];
        FRF = exp(1i * angle(Fopt * FBB'));
        mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];

        % measure elapsed time
        elapsed_time = toc(start_time);

        
        index = index + 1;
        cost(index) = real(mynorm(end));
        time(index) = elapsed_time;        
    end
    
    
    stats.cost = cost;
    stats.time = time;    
end