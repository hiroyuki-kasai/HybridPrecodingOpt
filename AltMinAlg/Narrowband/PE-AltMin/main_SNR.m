clear,clc
close all

Ns = 6;
file_str = sprintf('../../datasets/Ns=%d.mat', Ns);
load(file_str);

NRF = 8;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

for reali = 1:realization
    [ FRF, FBB ] = PE_AltMin( Fopt(:,:,reali), NRF);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = PE_AltMin( Wopt(:,:,reali), NRF);
    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
    
    fprintf('%04d: %s\n', reali, real(R(1,reali)));
end

plot(SNR_dB,sum(R,2)/realization,'g-->','LineWidth',1.5);
% plot(SNR_dB,sum(R,2)/realization,'Marker','>','LineWidth',1.5,'Color',[0 0.447058826684952 0.74117648601532]);
grid on
hold on
xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')