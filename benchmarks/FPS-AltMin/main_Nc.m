clear,clc

Nt = 144;
Nr = 16;
Ns = 4;
NRF = 4;

SNR_dB = 0;
SNR = 10.^(SNR_dB./10);
realization = 10000;

Nc_candi = 5:5:50;
CD = zeros(length(Nc_candi),1);

for i = 1:length(Nc_candi)
    Nc = Nc_candi(i)
    c = 1/sqrt(Nc)*exp(1i*[0:2*pi/Nc:2*pi-2*pi/Nc])';
    C = kron(eye(NRF),c);

    R = zeros(realization,1);
    parfor reali = 1:realization
        [ H,~,~,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns);
        
        [ FRF, FBB ] = my_AltMin_new( Fopt, C);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [ WRF, WBB ] = my_AltMin_new( Wopt, C);
        R(reali) = log2(det(eye(Ns) + SNR * pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
    end
    CD(i) = mean(R);
end

plot(Nc_candi,CD,'m-^','LineWidth',1.5);