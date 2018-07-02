% simulation code for fast hybrid precoding algorithm for OFDM.
% Created by H.Kasai on July 01, 2018
% Modified from the original codes in 
% https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems

clear
clc
close all;



%% set parameters
measure_snr_flag = false;
realization = 1000;

algorithms = {'Optimal', 'OMP', 'PE-AltMin', 'Proposed'};
alg_num = length(algorithms);


L = 128;    % # of multi-carriers
Nt = 144;   % # of transmit antennas
Nr = 36;    % # of receive antennas
%Nt = 16*16; % # of transmit antennas
%Nr = 8*8;   % # of receive antennas

Nc = 5;     % # of clusters
Nray = 10;  % # of rays in each cluster

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx

gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H


if measure_snr_flag
    Ns = 4; 
    NRF_array = 6;    
    SNR_dB = -15:5:10;
else
    Ns = 6;     
    NRF_array = 6:1:11;
    SNR_dB = 0;
end

SNR = 10.^(SNR_dB./10);
snr_num = length(SNR);  
nrf_num = length(NRF_array);

%% set params
clear params;
params.Ns = Ns;
params.Nc = Nc;
params.Nray = Nray;
params.Nt = Nt;
params.Nr = Nr;
params.angle_sigma = angle_sigma;  
params.gamma = gamma;    
params.sigma = sigma;    
params.realization = realization;       
params.L = L;  


%% assign and initialize cells
alg_elapsed_time = cell(alg_num, nrf_num);
info_cell        = cell(alg_num, realization);
speceff          = cell(alg_num, nrf_num);
speceff_sd       = cell(alg_num, nrf_num);

for j=1:alg_num
    for k=1:nrf_num
        alg_elapsed_time{j,k} = 0;
    end
end


%% perform algorithms
for nrf_idx = 1:length(NRF_array)
    NRF = NRF_array(nrf_idx);
    params.NRF = NRF;
    
    fprintf('##################### NRF = %d ##################### \n', NRF);
    
    RM = zeros(alg_num, snr_num, realization);
            
    for reali = 1:realization
        fprintf('# %d/%d \n', reali, realization);
        
        
        %% generate matrices for OFDM
        [Fopt, Wopt, H, At, Ar, FRF_enc, FRF_dec] = generate_OFDM_matrices(params);
            
        
        %% perform actual algorithms
        for alg_idx = 1 : alg_num

            alg_name = algorithms{alg_idx};
            
            start_time = tic();

            switch alg_name

                case 'Optimal'

                    % Do nothing

                case 'MO-AltMin'

                    [FRFM, FBBM, info_cell{alg_idx,reali}] = MO_AltMin(Fopt, NRF, FRF_enc);
                    for l = 1:L
                        FBBM(:,:,l) = sqrt(Ns) * FBBM(:,:,l) / norm(FRFM * FBBM(:,:,l),'fro');
                    end
                    [WRFM, WBBM, ~] = MO_AltMin( Wopt, NRF, FRF_dec );

                case 'Proposed'   

                    [FRFM, FBBM, info_cell{alg_idx,reali} ] = fast_hybrid_precoding_wrapper_OFDM( Fopt, NRF, FRF_enc);
                    for l = 1:L
                        FBBM(:,:,l) = sqrt(Ns) * FBBM(:,:,l) / norm(FRFM * FBBM(:,:,l),'fro');
                    end
                    [WRFM, WBBM, ~] = fast_hybrid_precoding_wrapper_OFDM( Wopt, NRF, FRF_dec);

                case 'PE-AltMin'  

                    [FRF, FBB, info_cell{alg_idx,reali}] = PE_AltMin_HK(Fopt, NRF, FRF_enc);
                    for l = 1:L
                        FBB(:,:,l) = sqrt(Ns) * FBB(:,:,l) / norm(FRF * FBB(:,:,l),'fro');
                    end
                    [WRF, WBB, ~] = PE_AltMin_HK( Wopt, NRF, FRF_dec ); 

                case 'OMP'

                    [FRFO, FBBO, info_cell{alg_idx,reali}] = OMP_HK(Fopt, NRF, At);
                    for l = 1:L
                        FBBO{l} = sqrt(Ns) * FBBO{l} / norm(FRFO * FBBO{l},'fro');
                    end
                    [WRFO, WBBO, ~] = OMP_HK(Wopt, NRF, Ar);                

            end
            
            elapsed_time = toc(start_time);
            alg_elapsed_time{alg_idx, nrf_idx} = alg_elapsed_time{alg_idx, nrf_idx} + elapsed_time;

 
            %% Calculate spectral efficiency
            for l = 1:L
                for s = 1:snr_num
                    if strcmp(alg_name, 'MO-AltMin') || strcmp(alg_name, 'Proposed')
                        RM(alg_idx,s,reali) = RM(alg_idx,s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFM * WBBM(:,:,l)) * H(:,:,l) * FRFM * FBBM(:,:,l) * FBBM(:,:,l)' * FRFM' * H(:,:,l)' * WRFM * WBBM(:,:,l)))/L;
                    elseif strcmp(alg_name, 'PE-AltMin') 
                        RM(alg_idx,s,reali) = RM(alg_idx,s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB(:,:,l)) * H(:,:,l) * FRF * FBB(:,:,l) * FBB(:,:,l)' * FRF' * H(:,:,l)' * WRF * WBB(:,:,l)))/L;
                    elseif strcmp(alg_name, 'OMP')
                        RM(alg_idx,s,reali) = RM(alg_idx,s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFO * WBBO{l}) * H(:,:,l) * FRFO * FBBO{l} * FBBO{l}' * FRFO' * H(:,:,l)' * WRFO * WBBO{l}))/L;
                    elseif strcmp(alg_name, 'Optimal')
                        RM(alg_idx,s,reali) = RM(alg_idx,s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,l)) * H(:,:,l) * Fopt(:,:,l) * Fopt(:,:,l)' * H(:,:,l)' * Wopt(:,:,l)))/L;                
                    end
                end
            end
            

            fprintf('[%s] \t%04d: %5.2f\n', alg_name, reali, real(RM(alg_idx,snr_num,reali)));
        
        end % end of algorithms
    end % end of realizations

        
    fprintf('\n');
    for alg_idx = 1 : alg_num
        alg_name = algorithms{alg_idx};
        speceff{alg_idx, nrf_idx}  = sum(real(RM(alg_idx,:,:)),3)/realization;
        speceff_sd{alg_idx, nrf_idx} = std(real(RM(alg_idx,:,:)),0,3);
        
        snr = speceff{alg_idx, nrf_idx};
        
        for s = 1:snr_num
            fprintf('## [%s] %d[dB]: \t%5.2f, %5.2f [sec]\n', alg_name, SNR_dB(s), snr(s), alg_elapsed_time{alg_idx, nrf_idx});
        end

    end
end 



%% plotting
fs = 20;
linewidth = 2;

style = {'-.+','-d','-*','-o'};
line_color = {[255, 128, 0], [76, 153, 0], [0, 0, 255], [255, 0, 0]};  

% generate legend string
legend_str = cell(length(algorithms),1);
for alg_idx = 1 : length(algorithms)
    legend_str{alg_idx} = algorithms{alg_idx};
end

if measure_snr_flag
    if snr_num == 1
        return;
    end
    
    % Spectral efficiency for each NRF
    for nrf_idx = 1:length(NRF_array)
        NRF = NRF_array(nrf_idx);
        
        title_str = sprintf('Spectral Efficiency for SNR (NRF=%d)', NRF);

        figure 
        title(title_str)
        for alg_idx = 1 : alg_num
            se = zeros(1,nrf_num);

            snr = speceff{alg_idx, nrf_idx};

            plot(SNR_dB, snr, style{alg_idx}, 'LineWidth', linewidth, 'Color', line_color{alg_idx}/255)
            hold on
        end
        
        grid on
        hold off
        ax1 = gca;
        set(ax1,'FontSize',fs);
        xlabel('SNR [dB]')
        ylabel('Spectral efficiency (bits/s/Hz)')
        legend(legend_str)
        title(title_str)
    end
    
else
    if nrf_num == 1
        return;
    end
    
    % Spectral efficiency when NRF changes
    for s=1:snr_num
        title_str = sprintf('Spectral Efficiency when NRF changes (%d[dB])', SNR_dB(s));

        figure 
        title(title_str)
        for alg_idx = 1 : alg_num
            se = zeros(1,nrf_num);
            for nrf_idx = 1:nrf_num
                snr = speceff{alg_idx, nrf_idx};
                snr_for_s(nrf_idx) = snr(s);

                %snr_sd = speceff_sd{alg_idx, nrf_idx};
                %snrsd_for_s(nrf_idx) = snr_sd(s);            
            end
            plot(NRF_array, snr_for_s, style{alg_idx}, 'LineWidth', linewidth, 'Color', line_color{alg_idx}/255)
            hold on
        end
        grid on
        hold off
        ax1 = gca;
        set(ax1,'FontSize',fs);
        xlabel('SNR [dB]')
        ylabel('Spectral efficiency (bits/s/Hz)')
        legend(legend_str)
        title(title_str)
    end
end




