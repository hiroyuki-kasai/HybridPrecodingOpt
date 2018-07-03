function [ FRF,FBB ] = my_AltMin_new( Fopt, C )
[Nt, Ns] = size(Fopt);
NRF = size(C,2);
Nc = size(C,1)/NRF;
[~,~,V] = svd(Fopt);
FBB = [V';zeros(NRF-Ns,Ns)];

mynorm = [Inf,0];
while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
    [alpha, value, s] = alpha_opt_new(real(Fopt*FBB'*C'));
    S = reshape(s,[Nt,Nc*NRF]);
    mynorm(1) = value + norm(imag(Fopt*FBB'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Fopt'*S*C,NRF);
    FBB = V*U';
    mynorm(2) = norm(Fopt*FBB'*C' - alpha*S,'fro')^2;
end

FBB = alpha*FBB;
FRF = S*C;
end