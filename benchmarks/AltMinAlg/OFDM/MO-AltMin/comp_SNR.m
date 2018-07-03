% This code is specialized for MO-AltMin, using parallel in MO-AltMin.m to speed up
clear
%clc
close all;


algorithms = {'MO_AltMin', 'MO_AltMin_new'};


Ns = 4; % # of streams
NRF = 5;

Nc = 5; % # of clusters
Nray = 10; % # of rays in each cluster

Nt = 144; % # of transmit antennas
Nr = 36; % # of receive antennas

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx

gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H

realization = 20;
count = 0;
K = 128;
%K = 12;

SNR_dB = -15:5:10;
SNR = 10.^(SNR_dB./10);
smax = length(SNR);

%% Init
Fopt_cell    = cell(realization,1);
Wopt_cell    = cell(realization,1);
H_cell       = cell(realization,1);
FRF_enc_cell = cell(realization,1);
FRF_dec_cell = cell(realization,1);



for reali = 1:realization
      
    Ht = zeros(Nr, Nt, Nc);
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2);
        AoA_m = unifrnd(0,2*pi,1,2);
        
        AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
        AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
        AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
        
        %Ht(:,:,c) = zeros(Nr,Nt);
        for j = 1:Nray
            temp = (c-1)*Nray+j;
            At(:,temp) = array_response(AoD(1,j),AoD(2,j),Nt);
            Ar(:,temp) = array_response(AoA(1,j),AoA(2,j),Nr);
            alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
            Ht(:,:,c) = Ht(:,:,c) + alpha * Ar(:,temp) * At(:,temp)';
        end
    end
    
    Fopt = zeros(Nt, Ns, K);
    Wopt = zeros(Nr, Ns, K);
    H = zeros(Nr, Nt, K);
    for k = 1:K
        %H(:,:,k) = zeros(Nr,Nt);
        for c = 1:Nc
            H(:,:,k) = H(:,:,k) + Ht(:,:,c) * exp(-1i*2*pi/K*(k-1)*(c-1));
        end
        H(:,:,k) = H(:,:,k) * gamma;
        if(rank(H(:,:,k))>=Ns)
            count = count + 1;
            
            [U,S,V] = svd(H(:,:,k));
            Fopt(:,:,k) = V([1:Nt],[1:Ns]);
            Wopt(:,:,k) = U([1:Nr],[1:Ns]);
        end
    end
    
    Fopt_cell{reali} = Fopt;
    Wopt_cell{reali} = Wopt;
    H_cell{reali} = H;    
    

    Nt_enc = size(Fopt, 1);
    FRF_enc_cell{reali} = exp( 1i*unifrnd(0,2*pi, Nt_enc, NRF) );  
    
    Nt_dec = size(Wopt, 1);
    FRF_dec_cell{reali} = exp( 1i*unifrnd(0,2*pi, Nt_dec, NRF) );     
end

figure
fs = 20;
linewidth = 2;
if 0

%% Conventional
RM = zeros(smax,realization);
tic
for reali = 1:realization
    
    
    %% MO-AltMin
    [ FRFM, FBBM, info_enc ] = MO_AltMin( Fopt_cell{reali}, NRF, FRF_enc_cell{reali} );
    for k = 1:K
        FBBM(:,:,k) = sqrt(Ns) * FBBM(:,:,k) / norm(FRFM * FBBM(:,:,k),'fro');
    end
    [ WRFM, WBBM, ~ ] = MO_AltMin( Wopt_cell{reali}, NRF, FRF_dec_cell{reali} );
    
    H = H_cell{reali};

    %% Calculate the spectral efficiency
    for k = 1:K
        for s = 1:smax
            RM(s,reali) = RM(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFM * WBBM(:,:,k)) * H(:,:,k) * FRFM * FBBM(:,:,k) * FBBM(:,:,k)' * FRFM' * H(:,:,k)' * WRFM * WBBM(:,:,k)))/K;
        end
    end
    
    fprintf('%04d: %s\n', reali, real(RM(1,reali)));
end
toc

sum(real(RM),2)/realization


plot(SNR_dB,sum(real(RM),2)'/realization,'-b','LineWidth',linewidth)
hold on;
end




%% Proposed
RM = zeros(smax,realization);
tic
for reali = 1:realization
    
    [ FRFM, FBBM, info_enc_HK ] = MO_AltMin_OFDM_HK( Fopt_cell{reali}, NRF, FRF_enc_cell{reali} );
    for k = 1:K
        FBBM(:,:,k) = sqrt(Ns) * FBBM(:,:,k) / norm(FRFM * FBBM(:,:,k),'fro');
    end
    [ WRFM, WBBM, ~ ] = MO_AltMin_OFDM_HK( Wopt_cell{reali}, NRF, FRF_dec_cell{reali} );

    H = H_cell{reali};    
    
    %% Calculate the spectral efficiency
    for k = 1:K
        for s = 1:smax
            RM(s,reali) = RM(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFM * WBBM(:,:,k)) * H(:,:,k) * FRFM * FBBM(:,:,k) * FBBM(:,:,k)' * FRFM' * H(:,:,k)' * WRFM * WBBM(:,:,k)))/K;
        end
    end
    
    fprintf('%04d: %s\n', reali, real(RM(1,reali)));
end
toc

sum(real(RM),2)/realization



plot(SNR_dB,sum(real(RM),2)'/realization,'-r','LineWidth',linewidth)
grid on
hold off
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('SNR [dB]')
ylabel('Spectral Efficiency (bits/s/Hz)')
legend('Conventional','Proposed')


%%
figure
plot(info_enc.cost,'-b','LineWidth',linewidth); hold on;
plot(info_enc_HK.cost,'-r','LineWidth',linewidth); hold off;
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('Iteration')
ylabel('Cost')
legend('Conventional','Proposed')


%%
figure
plot(info_enc.time, info_enc.cost,'-b','LineWidth',linewidth); hold on;
plot(info_enc_HK.time, info_enc_HK.cost,'-r','LineWidth',linewidth); hold off;
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('Time [sec]')
ylabel('Cost')
legend('Conventional','Proposed')


