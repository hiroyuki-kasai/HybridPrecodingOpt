function [ H,At,Ar,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns)

    %Ns = 4; % # of streams

    Nc = 5; % # of clusters
    Nray = 10; % # of rays in each cluster

    %Nt = 144; % # of transmit antennas
    %Nr = 36; % # of receive antennas

    angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx

    gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
    sigma = 1; %according to the normalization condition of the H

    %realization = 1000;
    count = 0;

    %for reali = 1:realization
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
            count = count + 1;

            [U,S,V] = svd(H(:,:));
            Fopt(:,:) = V([1:Nt],[1:Ns]);
            Wopt(:,:) = U([1:Nr],[1:Ns]);
        end
    %end
    
end