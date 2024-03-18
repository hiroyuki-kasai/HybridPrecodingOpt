% demonstratioin code for fast hybrid precoding algorithm for OFDM.

clear
clc
close all;


L = 128;    % # of multi-carriers
Ns = 4;     % # of streams
NRF = 5;    % # of RF

Nc = 5;     % # of clusters
Nray = 10;  % # of rays in each cluster

Nt = 144;   % # of transmit antennas
Nr = 36;    % # of receive antennas

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1;  %according to the normalization condition of the H

realization = 10;


SNR_dB = -15:5:10;
SNR = 10.^(SNR_dB./10);
smax = length(SNR);

% set params for generate_OFDM_matrices
clear params;
params.Ns = Ns;
params.NRF = NRF;
params.Nc = Nc;
params.Nray = Nray;
params.Nt = Nt;
params.Nr = Nr;
params.angle_sigma = angle_sigma;  
params.gamma = gamma;    
params.sigma = sigma;    
params.realization = realization;       
params.L = L;  

RM = zeros(smax, realization);


%% main loop with realization times
tic
for reali = 1:realization
 

    %% generate matrices for OFDM
    [Fopt, Wopt, H, At, Ar, FRF_enc, FRF_dec] = generate_OFDM_matrices(params);
    
    
    %% perform fast hybrid precoding algorithm
    [FRFM, FBBM] = fast_hybrid_precoding_OFDM_wrapper( Fopt, NRF, 0 );
    for l = 1:L
        FBBM(:,:,l) = sqrt(Ns) * FBBM(:,:,l) / norm(FRFM * FBBM(:,:,l),'fro');
    end
    [WRFM, WBBM] = fast_hybrid_precoding_OFDM_wrapper( Wopt, NRF, 0 );

    
    %% Calculate the spectral efficiency
    for l = 1:L
        for s = 1:smax
            RM(s,reali) = RM(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFM * WBBM(:,:,l)) * H(:,:,l) * FRFM * FBBM(:,:,l) * FBBM(:,:,l)' * FRFM' * H(:,:,l)' * WRFM * WBBM(:,:,l)))/L;
        end
    end
    
    fprintf('%04d: %s\n', reali, real(RM(1,reali)));
end
toc

sum(RM,2)/realization


%% plotting
figure
fs = 20;
linewidth = 2;

plot(SNR_dB,sum(RM,2)/realization, '-or', 'LineWidth', linewidth)
grid on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('Fast hybrid precoding algorithm', 'Location', 'northwest')
