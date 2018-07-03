clear,clc

%Nt = 144;
%Nr = 16;
%Ns = 4;
%NRF = 4;
Nt = 144;
Nr = 36;
Ns = 8;
NRF = 12;

SNR_dB = -15:5:10;
SNR = 10.^(SNR_dB./10);
realization = 200;
smax = length(SNR); % enable the parallel

Nc = 30;
c = 1/sqrt(Nc)*exp(1i*[0:2*pi/Nc:2*pi-2*pi/Nc])';
C = kron(eye(NRF),c);

%parfor reali = 1:realization
for reali = 1:realization
    reali
    [ H,At,Ar,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns);
    
    [ FRFA, FBBA ] = AE_AltMin( Fopt, NRF);
    FBBA = sqrt(Ns) * FBBA / norm(FRFA * FBBA,'fro');
    [ WRFA, WBBA ] = AE_AltMin( Wopt, NRF);
    
    [ FRFO, FBBO ] = OMP( Fopt, NRF, At);
    FBBO = sqrt(Ns) * FBBO / norm(FRFO * FBBO,'fro');
    [ WRFO, WBBO ] = OMP( Wopt, NRF, Ar);

    [ FRF, FBB ] = my_AltMin_new( Fopt, C);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = my_AltMin_new( Wopt, C);

    for s = 1:smax
        RA(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRFA*WBBA) * H * FRFA * FBBA * FBBA' * FRFA' * H' * WRFA*WBBA));
        ROM(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRFO*WBBO) * H * FRFO * FBBO * FBBO' * FRFO' * H' * WRFO*WBBO));
        R(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
        RO(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(Wopt) * H * Fopt * Fopt' * H' * Wopt));
    end
end
plot(SNR_dB,sum(RO,2)/realization,'r-o','LineWidth',1.5);
hold on
grid on
plot(SNR_dB,sum(R,2)/realization,'m-^','LineWidth',1.5);
plot(SNR_dB,sum(RA,2)/realization,'Marker','>','Color',[0 0.498039215803146 0],'LineWidth',1.5);
plot(SNR_dB,sum(ROM,2)/realization,'b-v','LineWidth',1.5);
legend('Fully digital','Proposed FPS-AltMin','MO-AltMin [3]','OMP [2]')