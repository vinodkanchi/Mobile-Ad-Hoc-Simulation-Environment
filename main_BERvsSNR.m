% -------------------------------------------------------------------------
%% ##### Initial Clean-Up #####
% -------------------------------------------------------------------------
clear ; close all; clc;

% -------------------------------------------------------------------------
%% ##### Simulation Settings #####
% -------------------------------------------------------------------------
c_fading        = true;     % Enable Rayleigh Fading
c_fadeAllSC     = false;    % Enable Subcarrier Fading
c_normPwr       = false;    % Normalize Transmit Power
c_filt          = false;    % Enable FIR filtering to suppress noise at RX

% -------------------------------------------------------------------------
%% ##### SNR Variation #####
% -------------------------------------------------------------------------
EbN0_vec            = 10 : 3 : +25;
N_packetsPerEbN0    = 10;
N_ANT               = 8;
%
N_runs              = N_ANT*length(EbN0_vec)*N_packetsPerEbN0;
antRange            = N_packetsPerEbN0*length(EbN0_vec);
%
% Console
N_printStatus       = ceil(0.05 * N_runs);

% -------------------------------------------------------------------------
%% ##### Allocate Memory #####
% -------------------------------------------------------------------------
SNR_tmp     = zeros(1,N_runs);
SNR_tmp_ch  = zeros(1,N_runs);
Ps_ch_vec   = zeros(1,N_runs);
Pn_ch_vec   = zeros(1,N_runs);
Ps_vec      = zeros(1,N_runs);
Pn_vec      = zeros(1,N_runs);
BER_tmp     = zeros(1,N_runs);

% -------------------------------------------------------------------------
%% ##### Map Data #####
% -------------------------------------------------------------------------
[bits_tx1,bits_tx2,bits_tx3,bits_tx4,bits_tx5,bits_tx6,bits_tx7,bits_tx8] = f_generate2packets();
data1 = qammod(bits_tx1,4,'InputType','bit','UnitAveragePower',true);
data2 = qammod(bits_tx2,4,'InputType','bit','UnitAveragePower',true);
data3 = qammod(bits_tx3,4,'InputType','bit','UnitAveragePower',true);
data4 = qammod(bits_tx4,4,'InputType','bit','UnitAveragePower',true);
data5 = qammod(bits_tx5,4,'InputType','bit','UnitAveragePower',true);
data6 = qammod(bits_tx6,4,'InputType','bit','UnitAveragePower',true);
data7 = qammod(bits_tx7,4,'InputType','bit','UnitAveragePower',true);
data8 = qammod(bits_tx8,4,'InputType','bit','UnitAveragePower',true);

% -------------------------------------------------------------------------
temp=0;
%% ##### Perform Simulation #####
% -------------------------------------------------------------------------
fprintf("##### Simulation started. #####\n")
t_start = tic;

currentRun = 0;

% LOOP: Antenna
 for k = 1 : 1 : N_ANT

    % Set Active Antennas
    A_tx1 = (k > 0) * 1;
    A_tx2 = (k > 1) * 1;
    A_tx3 = (k > 2) * 1;
    A_tx4 = (k > 3) * 1;
    A_tx5 = (k > 4) * 1;
    A_tx6 = (k > 5) * 1;
    A_tx7 = (k > 6) * 1;
    A_tx8 = (k > 7) * 1;
    
    %
    if c_normPwr
        N_TX = k;
    else
        N_TX = 1;
    end
    
    
    % LOOP: EbN0
    for j = 1 : 1 : length(EbN0_vec)
        % Calculate Noise Power
        EbN0_lin = 10^(EbN0_vec(j)/10);
        N0 = 1 / EbN0_lin; % Eb = 1;
        sigma_n = sqrt(N0);
        
        % LOOP: Transmission
        for i = 1 : 1 : N_packetsPerEbN0
            %% ##### Count Runs #####
            currentRun = currentRun + 1;
            
            %% ##### Generate OFDM Symbols and corresponding Bursts #####
            % TX 1
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H1_ray] = ...
                f_SIM_OFDM_tx1(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx1 = f_SIM_burst_tx1(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx1 = A_tx1 * u_tx1;
            % TX 2
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H2_ray] = ...
                f_SIM_OFDM_tx2(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx2 = f_SIM_burst_tx2(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx2 = A_tx2 * u_tx2;
            % TX 3
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H3_ray] = ...
                f_SIM_OFDM_tx3(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx3 = f_SIM_burst_tx3(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx3 = A_tx3 * u_tx3;
            % TX 4
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H4_ray] = ...
                f_SIM_OFDM_tx4(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx4 = f_SIM_burst_tx4(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx4 = A_tx4 * u_tx4;
            
            % TX 5
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H5_ray] = ...
                f_SIM_OFDM_tx5(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx5 = f_SIM_burst_tx5(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx5 = A_tx5 * u_tx5;
            % TX 6
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H6_ray] = ...
                f_SIM_OFDM_tx6(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx6 = f_SIM_burst_tx6(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx6 = A_tx6 * u_tx6;
            % TX 7
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H7_ray] = ...
                f_SIM_OFDM_tx7(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx7 = f_SIM_burst_tx7(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx7 = A_tx7 * u_tx7;
            % TX 8
            [ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H8_ray] = ...
                f_SIM_OFDM_tx8(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,false);
            u_tx8 = f_SIM_burst_tx8(ts1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,N_TX,EbN0_lin);
            u_tx8 = A_tx8 * u_tx8;
            %
            clear ts1 y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16

            
            %% ##### Channel #####
            % Super-Impose Transmit Bursts
            u_tx = u_tx1 + u_tx2 + u_tx3 + u_tx4 + u_tx5 + u_tx6 + u_tx7 + u_tx8;
            % Add AWGN (Delay == 1)
            L_rx = 6e4;
            rx = complex(zeros(L_rx, 1));
            rx(1:length(u_tx)) = u_tx;
            u_n_ch = sigma_n/sqrt(2)*(randn(size(rx))+1j*randn(size(rx)));
            rx = rx + u_n_ch;
            % Calc SNR on Channel
            Pn_ch = mean(abs(u_n_ch).^2);
            Ps_ch = mean(abs(u_tx(8*192*20+1:end)).^2);                     %% 4 -> 8
%             Pn_ch_vec(currentRun) = Pn_ch; % Debug Only
%             Ps_ch_vec(currentRun) = Ps_ch; % Debug Only
            %
            SNR_tmp_ch(currentRun) = 10*log10(Ps_ch/Pn_ch);

            
            %% ##### Sync #####
            H1 = A_tx1 * sum(abs(H1_ray.^2));
            H2 = A_tx2 * sum(abs(H2_ray.^2));
            H3 = A_tx3 * sum(abs(H3_ray.^2));
            H4 = A_tx4 * sum(abs(H4_ray.^2));
            H5 = A_tx5 * sum(abs(H5_ray.^2));
            H6 = A_tx6 * sum(abs(H6_ray.^2));
            H7 = A_tx7 * sum(abs(H7_ray.^2));
            H8 = A_tx8 * sum(abs(H8_ray.^2));
            %
            if (H1 > H2) && (H1 > H3) && (H1 > H4) && (H1 > H5) && (H1 > H6) && (H1 > H7) && (H1 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx1(rx,0.9,20,c_filt,true,1);
            elseif (H2 > H1) && (H2 > H3) && (H2 > H4) && (H2 > H5) && (H2 > H6) && (H2 > H7) && (H2 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx2(rx,0.9,20,c_filt,true,1);
            elseif (H3 > H1) && (H3 > H2) && (H3 > H4) && (H3 > H5) && (H3 > H6) && (H3 > H7) && (H3 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx3(rx,0.9,20,c_filt,true,1);
            elseif (H4 > H1) && (H4 > H2) && (H4 > H3) && (H4 > H5) && (H4 > H6) && (H4 > H7) && (H4 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx4(rx,0.9,20,c_filt,true,1);
            elseif (H5 > H1) && (H5 > H2) && (H5 > H3) && (H5 > H4) && (H5 > H6) && (H5 > H7) && (H5 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx5(rx,0.9,20,c_filt,true,1);
            elseif (H6 > H1) && (H6 > H2) && (H6 > H3) && (H6 > H4) && (H6 > H5) && (H6 > H7) && (H6 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx6(rx,0.9,20,c_filt,true,1);
            elseif (H7 > H1) && (H7 > H2) && (H7 > H3) && (H7 > H4) && (H7 > H5) && (H7 > H6) && (H7 > H8)
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx7(rx,0.9,20,c_filt,true,1);
            else
                [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,Ps,Pn] = ...
                    f_SIM_sync_tx8(rx,0.9,20,c_filt,true,1);
            end
            
            % Calculate SNR after FIR
%             Ps_vec(currentRun) = Ps; % Debug only
%             Pn_vec(currentRun) = Pn; % Debug only
            SNR_tmp(currentRun) = 10*log10(Ps/Pn);
            %
            clear rx H1 H2 H3 H4 H5 H6 H7 H8 Ps Pn Ps_ch Pn_ch

            %% ##### Decode #####
            % OFDM
            d1 = f_SIM_OFDM_decode(d1);
            d2 = f_SIM_OFDM_decode(d2);
            d3 = f_SIM_OFDM_decode(d3);
            d4 = f_SIM_OFDM_decode(d4);
            d5 = f_SIM_OFDM_decode(d5);
            d6 = f_SIM_OFDM_decode(d6);
            d7 = f_SIM_OFDM_decode(d7);
            d8 = f_SIM_OFDM_decode(d8);
            d9 = f_SIM_OFDM_decode(d9);
            d10 = f_SIM_OFDM_decode(d10);
            d11 = f_SIM_OFDM_decode(d11);
            d12 = f_SIM_OFDM_decode(d12);
            d13 = f_SIM_OFDM_decode(d13);
            d14 = f_SIM_OFDM_decode(d14);
            d15 = f_SIM_OFDM_decode(d15);
            d16 = f_SIM_OFDM_decode(d16);           
            % OSTBC
            [sym1_rx,sym2_rx,sym3_rx,sym4_rx,sym5_rx,sym6_rx,sym7_rx,sym8_rx] = ...
                f_SIM_OSTBC_decode(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,...
                    H1_ray,H2_ray,H3_ray,H4_ray,H5_ray,H6_ray,H7_ray,H8_ray,...
                    A_tx1,A_tx2,A_tx3,A_tx4,A_tx5,A_tx6,A_tx7,A_tx8);
            % AGC
            sym1_rx = f_SIM_AGC(sym1_rx,true);
            sym2_rx = f_SIM_AGC(sym2_rx,true);
            sym3_rx = f_SIM_AGC(sym3_rx,true);
            sym4_rx = f_SIM_AGC(sym4_rx,true);
            sym5_rx = f_SIM_AGC(sym5_rx,true);
            sym6_rx = f_SIM_AGC(sym6_rx,true);
            sym7_rx = f_SIM_AGC(sym7_rx,true);
            sym8_rx = f_SIM_AGC(sym8_rx,true);
            % 4-QAM
            bits_rx1 = qamdemod(sym1_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx2 = qamdemod(sym2_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx3 = qamdemod(sym3_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx4 = qamdemod(sym4_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx5 = qamdemod(sym5_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx6 = qamdemod(sym6_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx7 = qamdemod(sym7_rx,4,'OutputType','bit','UnitAveragePower',true);
            bits_rx8 = qamdemod(sym8_rx,4,'OutputType','bit','UnitAveragePower',true);            
            %
            clear d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16
            clear sym1_rx sym2_rx sym3_rx sym4_rx sym5_rx sym6_rx sym7_rx sym8_rx

            %% ##### BER #####
            BER_1 = sum(abs(bits_tx1-bits_rx1)) / length(bits_tx1);
            BER_2 = sum(abs(bits_tx2-bits_rx2)) / length(bits_tx2);
            BER_3 = sum(abs(bits_tx3-bits_rx3)) / length(bits_tx3);
            BER_4 = sum(abs(bits_tx4-bits_rx4)) / length(bits_tx4);
            BER_5 = sum(abs(bits_tx5-bits_rx5)) / length(bits_tx5);
            BER_6 = sum(abs(bits_tx6-bits_rx6)) / length(bits_tx6);
            BER_7 = sum(abs(bits_tx7-bits_rx7)) / length(bits_tx7);
            BER_8 = sum(abs(bits_tx8-bits_rx8)) / length(bits_tx8);
            BER_tmp(currentRun) = 1/8 * (BER_1 + BER_2 + BER_3 + BER_4 + BER_5 + BER_6 + BER_7 + BER_8);

            % -----------------------------------------------------------------
            %% Print Status
            if (rem(currentRun,N_printStatus) == 0)
                f_printSimStat(currentRun,N_runs,EbN0_vec(j),...
                    BER_tmp(currentRun));
            end
        end % Transmission

    end % EbN0
%     
end % Antenna


% -------------------------------------------------------------------------
%% ##### Convert Elapsed Time #####
% -------------------------------------------------------------------------
t_elapsed = round(toc(t_start),3);
t_text    = f_convertElapsedTime(t_elapsed);
fprintf("\n-> Simulation completed in " + t_text + ".\n\n");


% -------------------------------------------------------------------------
%% ##### Save Results #####
% -------------------------------------------------------------------------
% Allocate Memory
BER = zeros(N_ANT,antRange);
SNR = zeros(N_ANT,antRange);
SNR_ch = zeros(N_ANT,antRange);
% Split Up Vectors
if N_ANT > 0
    BER(1,:) = BER_tmp(0*antRange+1 : 1*antRange);
    SNR(1,:) = SNR_tmp(0*antRange+1 : 1*antRange);
    SNR_ch(1,:) = SNR_tmp_ch(0*antRange+1 : 1*antRange);
end
if N_ANT > 1
    BER(2,:) = BER_tmp(1*antRange+1 : 2*antRange);
    SNR(2,:) = SNR_tmp(1*antRange+1 : 2*antRange);
    SNR_ch(2,:) = SNR_tmp_ch(1*antRange+1 : 2*antRange);
end
if N_ANT > 2
    BER(3,:) = BER_tmp(2*antRange+1 : 3*antRange);
    SNR(3,:) = SNR_tmp(2*antRange+1 : 3*antRange);
    SNR_ch(3,:) = SNR_tmp_ch(2*antRange+1 : 3*antRange);
end
if N_ANT > 3
    BER(4,:) = BER_tmp(3*antRange+1 : 4*antRange);
    SNR(4,:) = SNR_tmp(3*antRange+1 : 4*antRange);
    SNR_ch(4,:) = SNR_tmp_ch(3*antRange+1 : 4*antRange);
end
if N_ANT > 4
    BER(5,:) = BER_tmp(4*antRange+1 : 5*antRange);
    SNR(5,:) = SNR_tmp(4*antRange+1 : 5*antRange);
    SNR_ch(5,:) = SNR_tmp_ch(4*antRange+1 : 5*antRange);
end
if N_ANT > 5
    BER(6,:) = BER_tmp(5*antRange+1 : 6*antRange);
    SNR(6,:) = SNR_tmp(5*antRange+1 : 6*antRange);
    SNR_ch(6,:) = SNR_tmp_ch(5*antRange+1 : 6*antRange);
end
if N_ANT > 6
    BER(7,:) = BER_tmp(6*antRange+1 : 7*antRange);
    SNR(7,:) = SNR_tmp(6*antRange+1 : 7*antRange);
    SNR_ch(7,:) = SNR_tmp_ch(6*antRange+1 : 7*antRange);
end
if N_ANT > 7
    BER(8,:) = BER_tmp(7*antRange+1 : 8*antRange);
    SNR(8,:) = SNR_tmp(7*antRange+1 : 8*antRange);
    SNR_ch(8,:) = SNR_tmp_ch(7*antRange+1 : 8*antRange);
end
% Save Results
save("./SIM_results/resFile_BER.mat",'BER');
save("./SIM_results/resFile_SNR.mat",'SNR');
save("./SIM_results/resFile_SNR_ch.mat",'SNR_ch');
save("./SIM_results/resFile_EbN0.mat",'EbN0_vec');
save("./SIM_results/resFile_Pn_ch.mat",'Pn_ch_vec');  % Debug Only
save("./SIM_results/resFile_Ps_ch.mat",'Ps_ch_vec');  % Debug Only
save("./SIM_results/resFile_Pn.mat",'Pn_vec');        % Debug Only
save("./SIM_results/resFile_Ps.mat",'Ps_vec');        % Debug Only
%
fprintf("-> Saved results to ./SIM_results/...mat.\n\n")


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% ##### Plot Results #####
% -------------------------------------------------------------------------
%% Plot Settings
usr_marker      = '-x';
usr_LineWidth   = 2;
usr_MarkerSize  = 4;
usr_FontSize    = 10.5;

%% Calc Mean Vectors
EbN0_steps = length(EbN0_vec);
SNR_thres = -79;
if N_ANT > 0
    BER_1       = f_calc_meanBERvec(BER(1,:),N_packetsPerEbN0,EbN0_steps);
    SNR_1       = f_calc_meanSNRvec(SNR(1,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_1_ch    = f_calc_meanSNRvec(SNR_ch(1,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_1        = Ps_vec(0*N_runs/8+1 : 1*N_runs/8);
    Ps_ch_1     = Ps_ch_vec(0*N_runs/8+1 : 1*N_runs/8);
    Pn_1        = Pn_vec(0*N_runs/8+1 : 1*N_runs/8);
    Pn_ch_1     = Pn_ch_vec(0*N_runs/8+1 : 1*N_runs/8);
end
if N_ANT > 1
    BER_2       = f_calc_meanBERvec(BER(2,:),N_packetsPerEbN0,EbN0_steps);
    SNR_2       = f_calc_meanSNRvec(SNR(2,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_2_ch    = f_calc_meanSNRvec(SNR_ch(2,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_2        = Ps_vec(1*N_runs/8+1 : 2*N_runs/8);
    Ps_ch_2     = Ps_ch_vec(1*N_runs/8+1 : 2*N_runs/8);
    Pn_2        = Pn_vec(1*N_runs/8+1 : 2*N_runs/8);
    Pn_ch_2     = Pn_ch_vec(1*N_runs/8+1 : 2*N_runs/8);
end
if N_ANT > 2
    BER_3       = f_calc_meanBERvec(BER(3,:),N_packetsPerEbN0,EbN0_steps);
    SNR_3       = f_calc_meanSNRvec(SNR(3,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_3_ch    = f_calc_meanSNRvec(SNR_ch(3,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_3        = Ps_vec(2*N_runs/8+1 : 3*N_runs/8);
    Ps_ch_3     = Ps_ch_vec(2*N_runs/8+1 : 3*N_runs/8);
    Pn_3        = Pn_vec(2*N_runs/8+1 : 3*N_runs/8);
    Pn_ch_3     = Pn_ch_vec(2*N_runs/8+1 : 3*N_runs/8);
end
if N_ANT > 3
    BER_4       = f_calc_meanBERvec(BER(4,:),N_packetsPerEbN0,EbN0_steps);
    SNR_4       = f_calc_meanSNRvec(SNR(4,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_4_ch    = f_calc_meanSNRvec(SNR_ch(4,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_4        = Ps_vec(3*N_runs/8+1 : 4*N_runs/8);
    Ps_ch_4     = Ps_ch_vec(3*N_runs/8+1 : 4*N_runs/8);
    Pn_4        = Pn_vec(3*N_runs/8+1 : 4*N_runs/8);
    Pn_ch_4     = Pn_ch_vec(3*N_runs/8+1 : 4*N_runs/8);
end
if N_ANT > 4
    BER_5       = f_calc_meanBERvec(BER(5,:),N_packetsPerEbN0,EbN0_steps);
    SNR_5       = f_calc_meanSNRvec(SNR(5,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_5_ch    = f_calc_meanSNRvec(SNR_ch(5,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_5        = Ps_vec(4*N_runs/8+1 : 5*N_runs/8);
    Ps_ch_5     = Ps_ch_vec(4*N_runs/8+1 : 5*N_runs/8);
    Pn_5        = Pn_vec(4*N_runs/8+1 : 5*N_runs/8);
    Pn_ch_5     = Pn_ch_vec(4*N_runs/8+1 : 5*N_runs/8);
end
if N_ANT > 5
    BER_6       = f_calc_meanBERvec(BER(6,:),N_packetsPerEbN0,EbN0_steps);
    SNR_6       = f_calc_meanSNRvec(SNR(6,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_6_ch    = f_calc_meanSNRvec(SNR_ch(6,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_6        = Ps_vec(5*N_runs/8+1 : 6*N_runs/8);
    Ps_ch_6     = Ps_ch_vec(5*N_runs/8+1 : 6*N_runs/8);
    Pn_6        = Pn_vec(5*N_runs/8+1 : 6*N_runs/8);
    Pn_ch_6     = Pn_ch_vec(5*N_runs/8+1 : 6*N_runs/8);
end
if N_ANT > 6
    BER_7       = f_calc_meanBERvec(BER(7,:),N_packetsPerEbN0,EbN0_steps);
    SNR_7       = f_calc_meanSNRvec(SNR(7,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_7_ch    = f_calc_meanSNRvec(SNR_ch(7,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_7        = Ps_vec(6*N_runs/8+1 : 7*N_runs/8);
    Ps_ch_7     = Ps_ch_vec(6*N_runs/8+1 : 7*N_runs/8);
    Pn_7        = Pn_vec(6*N_runs/8+1 : 7*N_runs/8);
    Pn_ch_7     = Pn_ch_vec(6*N_runs/8+1 : 7*N_runs/8);
end
if N_ANT > 7
    BER_8       = f_calc_meanBERvec(BER(8,:),N_packetsPerEbN0,EbN0_steps);
    SNR_8       = f_calc_meanSNRvec(SNR(8,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    SNR_8_ch    = f_calc_meanSNRvec(SNR_ch(8,:),N_packetsPerEbN0,EbN0_steps,SNR_thres);
    Ps_8        = Ps_vec(7*N_runs/8+1 : 8*N_runs/8);
    Ps_ch_8     = Ps_ch_vec(7*N_runs/8+1 : 8*N_runs/8);
    Pn_8        = Pn_vec(7*N_runs/8+1 : 8*N_runs/8);
    Pn_ch_8     = Pn_ch_vec(7*N_runs/8+1 : 8*N_runs/8);
end

%% Plot: EbN0 vs. SNR
figure()
semilogy(EbN0_vec,EbN0_vec,usr_marker,'LineWidth',1.5*usr_LineWidth,'MarkerSize',1.5*usr_MarkerSize);
hold on;
if N_ANT > 0
   semilogy(EbN0_vec,SNR_1,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1"])
end
if N_ANT > 1
   semilogy(EbN0_vec,SNR_2,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1";"SNR_2"])
end
if N_ANT > 2
   semilogy(EbN0_vec,SNR_3,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3"])
end
if N_ANT > 3
   semilogy(EbN0_vec,SNR_4,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize); 
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3";"SNR_4"])
end
if N_ANT > 4
   semilogy(EbN0_vec,SNR_5,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3";"SNR_4";"SNR_5"])
end
if N_ANT > 5
   semilogy(EbN0_vec,SNR_6,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3";"SNR_4";"SNR_5";"SNR_6"])
end
if N_ANT > 6
   semilogy(EbN0_vec,SNR_7,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3";"SNR_4";"SNR_5";"SNR_6";"SNR_7"])
end
if N_ANT > 7
   semilogy(EbN0_vec,SNR_8,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize); hold off;
   legend(["EbN0";"SNR_1";"SNR_2";"SNR_3";"SNR_4";"SNR_5";"SNR_6";"SNR_7";"SNR_8"])
end
grid()
set(gca,'FontSize',usr_FontSize)
title("SNR vs. EbN0 per Transmission")
ylabel("SNR / EbN0")
xlabel("EbN0 in dB");

%% Plot: BER vs. SNR
figure()
hold on
if N_ANT > 0
   semilogy(EbN0_vec-13,BER_1,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize); 
   legend("1 TX")
end
if N_ANT > 1
   semilogy(EbN0_vec-13,BER_2,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX"])
end
if N_ANT > 2
   semilogy(EbN0_vec-13,BER_3,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX";"3 TX"])
end
if N_ANT > 3
   semilogy(EbN0_vec-13,BER_4,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX";"3 TX";"4 TX"])
end
if N_ANT > 4
   semilogy(EbN0_vec-13,BER_5,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX"])
end
if N_ANT > 5
   semilogy(EbN0_vec-13,BER_6,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX"])
end
if N_ANT > 6
   semilogy(EbN0_vec-13,BER_7,usr_marker,'LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX"])
end
if N_ANT > 7
   semilogy(EbN0_vec-13,BER_8,'-b','LineWidth',usr_LineWidth,'MarkerSize',usr_MarkerSize); 
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX";"8 TX"])
end
hold off;
grid()
set(gca,'FontSize',usr_FontSize)
title("BER vs. Transmit Power")
ylabel("BER")
xlabel("SNR in dB");


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% ##### Clean-Up #####
% -------------------------------------------------------------------------
clear A_tx1 A_tx2 A_tx3 A_tx4 A_tx5 A_tx6 A_tx7 A_tx8 N_TX N_ANT data_scale
clear u_tx1 u_tx2 u_tx3 u_tx4 u_tx5 u_tx6 u_tx7 u_tx8 u_tx u_n_ch rx L_rx
clear c_fading c_fadeAllSC c_filt c_normPwr H1_ray H2_ray H3_ray H4_ray H5_ray H6_ray H7_ray H8_ray 
clear bits_tx1 bits_tx2 bits_tx3 bits_tx4 bits_tx5 bits_tx6 bits_tx7 bits_tx8 data1 data2 data3 data4 data5 data6 data7 data8 N_data
clear bits_rx1 bits_rx2 bits_rx3 bits_rx4 bits_rx5 bits_rx6 bits_rx7 bits_rx8 BER_1 BER_2 BER_3 BER_4 BER_5 BER_6 BER_7 BER_8
clear i j k currentRun N_runs t_start t_elapsed t_h t_min t_sec
clear sigma_n N0 EbN0_lin N_packetsPerEbN0 EbN0_cntr ant_cntr antRange
clear BER_tmp SNR_tmp SNR_tmp_ch SNR_thres EbN0_steps
clear usr_FontSize usr_LineWidth usr_marker usr_MarkerSize