function [FRF, FBB, stats] = fast_hybrid_precoding_NB_wrapper(Fopt, NRF, FRF_in, options_in)
% Wapper function for fast hybrid_precoding Narrowband algorithm.
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


    if ~isfield(options_in, 'solver')
        options.solver = 'cg';
    else
        options.solver = options_in.solver;
    end  
    
    if ~isfield(options_in, 'func_tolerance')
        options.func_tolerance = 1e-1;
    else
        options.func_tolerance = options_in.func_tolerance;
    end         
    

    [Nt, Ns] = size(Fopt);
    
    % FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    if FRF_in == 0
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    else
        FRF = FRF_in;
    end 
    
    % calculate initial cost and FBB
    FBB = pinv(FRF) * Fopt;
    init_cost = norm(Fopt(:,:) - FRF * FBB(:,:),'fro')^2;

    %options.solver = 'cg'; 
    %options.func_tolerance = options_in.func_tolerance;      
    
    [FRF,FBB, info] = fast_hybrid_precoding_NB(Fopt, FRF, FBB, options);
    
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