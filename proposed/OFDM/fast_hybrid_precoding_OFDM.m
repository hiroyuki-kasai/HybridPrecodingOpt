function [FRF, FBB, info] = fast_hybrid_precoding_OFDM(Fopt, FRF_init, K, options)
% Fast hybrid precoding OFDM algorithm.
%
% Inputs:
%       Fopt        optimul fully digital precoder matrix
%       FRF_init    initial matrix for analog RF precoder
%       K           number of malti-carriers
%       options     some options

% Output:
%       FRF         calculated matrix for analog RF precoder
%       FBB         calculated matrix for digital baseband precoder
%       info        statistics
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


    if ~isfield(options, 'solver')
        options.solver = 'cg';
    end  
    
    if ~isfield(options, 'func_tolerance')
        options.func_tolerance = 1e-1;
    end  
    
    
    [Nt, NRF] = size(FRF_init);
    %K = size(FBB,3);

    % set manifold
    manifold = complexcirclefactory(Nt * NRF);
    problem.M = manifold;

    f_prev = 1000;

    problem.cost = @(x) mycost(x);
    function f_val = mycost(x)
        
        FRF_cur = reshape(x, Nt, NRF);
        invFRF_cur = pinv(FRF_cur);
        
        f_val = 0;
        for k = 1:K
            FBB(:,:,k) = invFRF_cur * Fopt(:,:,k);
            diff_mat = Fopt(:,:,k) - FRF_cur * FBB(:,:,k);
            f_val = f_val + diff_mat(:)' * diff_mat(:);
        end         
        
    end   
    
    problem.egrad = @(x) mygrad(x);
    function g = mygrad(x)
        
        FRF_cur = reshape(x, Nt ,NRF);
        invFRF_cur = pinv(FRF_cur);
        
        g_mat = 0;
        for k = 1:K
            FBB(:,:,k) = invFRF_cur * Fopt(:,:,k); 
            gg = (Fopt(:,:,k) - FRF_cur*FBB(:,:,k)) * FBB(:,:,k)';         
            g_mat = g_mat + gg;
        end  
        
        % perform vectorization
        g = -2 * g_mat(:);

    end  

    maopt_options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, x, info_stop, last)
        stopnow = 0;
        
        f_curr = mycost(x);
        if abs(f_curr - f_prev) < options.func_tolerance
            stopnow = 1;
        end
        f_prev = f_curr;
    end


    % checkgradient(problem);
    warning('off', 'manopt:getHessian:approx');

    maopt_options.verbosity = 0;
    
    if strcmp(options.solver, 'sd')    
        [x, cost, info, maopt_options] = steepestdescent(problem, FRF_init(:), maopt_options);
    elseif strcmp(options.solver, 'cg')
        [x, cost, info, maopt_options] = conjugategradient(problem, FRF_init(:), maopt_options);
    else
        [x, cost, info, maopt_options] = trustregions(problem, FRF_init(:), maopt_options);
    end
    
    % perform inverse-vectorization
    FRF = reshape(x, Nt, NRF);
end