% demonstratioin code for fast hybrid precoding algorithm for Narrowband.
% Created by H.Kasai on July 01, 2018
% Modified from the original codes in 
% https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems

close all;
clear
clc


% load('Ns=3.mat');
Ns = 3;
filestr = sprintf('../../benchmarks/AltMinAlg/datasets/Ns=%d.mat', Ns);
load(filestr);


NRF = 3;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

realization = 5;

tic
for reali = 1:realization
    
    %% perform fast hybrid precoding algorithm
    [FRF, FBB] = fast_hybrid_precoding_wrapper_NB( Fopt(:,:,reali), NRF, 0);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [WRF, WBB] = fast_hybrid_precoding_wrapper_NB( Wopt(:,:,reali), NRF, 0);

    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
    
    fprintf('%04d: %s\n', reali, real(R(1,reali)));
end
toc

sum(R,2)/realization


%% plotting
figure
fs = 20;
linewidth = 2;
plot(SNR_dB,sum(R,2)/realization, '-or', 'LineWidth', linewidth)
grid on
hold on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')

