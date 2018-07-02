function [y, FBB, cost] = sig_manif_NB_HK(Fopt, FRF, FBB)
    [Nt, NRF] = size(FRF);

    manifold = complexcirclefactory(Nt*NRF);
    %manifold = obliquecomplexfactory(1,Nt*NRF,true);
    problem.M = manifold;

    % problem.cost  = @(x) norm( Fopt - reshape(x,Nt,NRF) * FBB,'fro')^2;
    % problem.egrad = @(x) -2 * kron(conj(FBB), eye(Nt)) * (Fopt(:) - kron(FBB.', eye(Nt)) * x);
    f = Fopt(:);
    A = kron(FBB.', eye(Nt));

    %problem.cost  = @(x) (f-A*x)'*(f-A*x);
    %problem.egrad = @(x) -2*A'*(f-A*x);

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

        %FBB = pinv(FRF) * Fopt;
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

    % checkgradient(problem);
    warning('off', 'manopt:getHessian:approx');

    options.verbosity = 0;

    [x,cost,info,options] = conjugategradient(problem,FRF(:), options);
    % [x,cost,info,options] = trustregions(problem, FRF(:), options);
    % info.iter
    y = reshape(x,Nt,NRF);

end