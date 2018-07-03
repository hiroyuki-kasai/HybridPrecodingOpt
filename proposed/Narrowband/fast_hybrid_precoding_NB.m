function [y, FBB, info] = fast_hybrid_precoding_NB(Fopt, FRF_init, FBB, options)
% Fast hybrid precoding Narrowband algorithm.
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

    manifold = complexcirclefactory(Nt*NRF);
    %manifold = obliquecomplexfactory(1,Nt*NRF,true);
    problem.M = manifold;

    % problem.cost  = @(x) norm( Fopt - reshape(x,Nt,NRF) * FBB,'fro')^2;
    % problem.egrad = @(x) -2 * kron(conj(FBB), eye(Nt)) * (Fopt(:) - kron(FBB.', eye(Nt)) * x);
    f = Fopt(:);
    A = kron(FBB.', eye(Nt));

    %problem.cost  = @(x) (f-A*x)'*(f-A*x);
    %problem.egrad = @(x) -2*A'*(f-A*x);
    

    f_prev = 1000;    

    problem.cost = @mycost;
    function [f_value] = mycost(x)

         FRF_cur = reshape(x,Nt,NRF);
%         FBB = pinv(FRF_cur) * Fopt;
%         A = kron(FBB.', eye(Nt));
% 
%         f_value = (f-A*x)'*(f-A*x);
        
        diff = Fopt - FRF_cur * FBB;
        f_value = diff(:)' * diff(:);
    end

    problem.egrad = @myegrad;
    function [g] = myegrad(x)

        %FBB = pinv(FRF_init) * Fopt;
        FRF_cur = reshape(x,Nt,NRF);
        FBB = pinv(FRF_cur) * Fopt;

        if 0
            A = kron(FBB.', eye(Nt));
            g = -2*A'*(f-A*x);
        else
            %x_mat = reshape(x,Nt,NRF);
            gg = -2*(Fopt-FRF_cur*FBB)*FBB';
            g = reshape(gg, Nt*NRF, 1);
        end
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

    maopt_options.verbosity = 0;
    
    if strcmp(options.solver, 'sd')    
        [x, cost, info, maopt_options] = steepestdescent(problem, FRF_init(:), maopt_options);
    elseif strcmp(options.solver, 'cg')
        [x, cost, info, maopt_options] = conjugategradient(problem, FRF_init(:), maopt_options);
    else
        [x, cost, info, maopt_options] = trustregions(problem, FRF_init(:), maopt_options);
    end
    
    y = reshape(x, Nt, NRF);

end