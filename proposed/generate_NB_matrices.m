function [Fopt, Wopt, H, At, Ar, FRF_enc, FRF_dec] = generate_NB_matrices(params)

    Ns = params.Ns;
    NRF = params.NRF;
    Nc = params.Nc;
    Nray = params.Nray;
    Nt = params.Nt;
    Nr = params.Nr;
    angle_sigma = params.angle_sigma;  
    gamma = params.gamma;    
    sigma = params.sigma;    


%     Ht = zeros(Nr, Nt, Nc);
%     for c = 1:Nc
%         AoD_m = unifrnd(0,2*pi,1,2);
%         AoA_m = unifrnd(0,2*pi,1,2);
% 
%         AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
%         AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
%         AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
%         AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
% 
%         for j = 1:Nray
%             temp = (c-1)*Nray+j;
%             At(:,temp) = array_response(AoD(1,j),AoD(2,j),Nt);
%             Ar(:,temp) = array_response(AoA(1,j),AoA(2,j),Nr);
%             alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
%             Ht(:,:,c) = Ht(:,:,c) + alpha * Ar(:,temp) * At(:,temp)';
%         end
%     end
% 
%     Fopt = zeros(Nt, Ns);
%     Wopt = zeros(Nr, Ns);
%     H = zeros(Nr, Nt);
%     for c = 1:Nc
%         H(:,:) = H(:,:) + Ht(:,:,c) * exp(-1i*2*pi/(l-1)*(c-1));
%     end
%     H(:,:) = H(:,:) * gamma;
%     if(rank(H(:,:,l))>=Ns)
%         count = count + 1;
% 
%         [U,S,V] = svd(H(:,:,l));
%         Fopt(:,:) = V(1:Nt, 1:Ns);
%         Wopt(:,:) = U(1:Nr, 1:Ns);
%     end
    
    
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2);
        AoA_m = unifrnd(0,2*pi,1,2);
        
        AoD(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(1),angle_sigma);
        AoD(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(2),angle_sigma);
        AoA(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    end
    
    H(:,:) = zeros(Nr,Nt);
    for j = 1:Nc*Nray
        At(:,j) = array_response(AoD(1,j),AoD(2,j),Nt); %UPA array response
        Ar(:,j) = array_response(AoA(1,j),AoA(2,j),Nr);
        alpha(j) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
        H(:,:) = H(:,:) + alpha(j) * Ar(:,j) * At(:,j)';
    end
    H(:,:) = gamma * H(:,:);
    
    if(rank(H(:,:))>=Ns)
    
        [U,S,V] = svd(H(:,:));
        Fopt(:,:) = V(1:Nt, 1:Ns);
        Wopt(:,:) = U(1:Nr, 1:Ns);
    end    


    Nt_enc = size(Fopt, 1);
    FRF_enc = exp( 1i*unifrnd(0,2*pi, Nt_enc, NRF) );  

    Nt_dec = size(Wopt, 1);
    FRF_dec = exp( 1i*unifrnd(0,2*pi, Nt_dec, NRF) );     
end

