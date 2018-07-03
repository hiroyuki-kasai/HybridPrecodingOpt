function [FRF, FBB, stats] = fast_hybrid_precoding_wrapper_OFDM(Fopt, NRF, FRF_in)
% Wapper function for fast hybrid_precoding OFDM algorithm.
%
% Inputs:
%       Fopt        optimul fully digital precoder matrix
%       NRF         number of radio frequency (RF) 
%       FRF_in      initial matrix for analog RF precoder
% Output:
%       FRF         calculated matrix for analog RF precoder
%       FBB         calculated matrix for digital baseband precoder
%       stats       statistics
%
% Reference:
%       Hiroyuki Kasai, 
%       "Fast optimization algorithm on complex oblique manifold for
%       hybrid precoding in Millimeter Wave MIMO systems,"
%       arXiv, 2018.
%
%
% Created by H.Kasai on July 01, 2018
% Modified from the original codes in 
% https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems


    [Nt, Ns, K] = size(Fopt);

    if FRF_in == 0
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    else
        FRF = FRF_in;
    end  
    
    % calculate initial cost and FBB
    FBB = zeros(NRF, Ns, K);
    init_cost = 0;
    for k = 1:K      
        FBB(:,:,k) = pinv(FRF) * Fopt(:,:,k);
        init_cost = init_cost + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
    end

    options.solver = 'cg'; 
    options.func_tolerance = 1e-1;     
    [FRF, FBB, info] = fast_hybrid_precoding_OFDM(Fopt, FRF, K, options);
    
    % store stats infos
    len = length(info);
    cost = zeros(1, len+1);
    time = zeros(1, len+1);
    cost(1) = init_cost;
    time(1) = 0;    
    for j=1:len
        cost(1+j) = info(j).cost;
        time(1+j) = info(j).time;
    end
    
    stats.cost = cost;
    stats.time = time;
end