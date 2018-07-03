clear,
clc

load('../../datasets/Ns=3.mat');


Ns = 3;

NRF = 3;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

realization = 10;

tic
for reali = 1:realization
    
    [ FRF, FBB ] = MO_AltMin( Fopt(:,:,reali), NRF, 0);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = MO_AltMin( Wopt(:,:,reali), NRF, 0);

    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
    
    fprintf('%04d: %s\n', reali, real(R(1,reali)));
end
toc

plot(SNR_dB,sum(R,2)/realization,'k-p','LineWidth',1.5)
grid on
hold on
xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')

sum(R,2)/realization