% -------------------------------------------------------------------------
%% ##### Initial Clean-Up #####
% -------------------------------------------------------------------------
clear all; close all; clc;

% -------------------------------------------------------------------------
%% ##### SETTINGS #####
% -------------------------------------------------------------------------
% SIM Result Files
file_nbr    = 1;
% Simulation Settings
SNR_thres   = -79;
% Plot Settings
usr_mrkr    = '-x';     % Marker
usr_LW      = 2;        % Line-Width
usr_MS      = 4;        % Marker-Size
usr_FS      = 10.5;     % Font-Size
c_xlim      = false;
usr_xlim    = [7; 25]; % x-axis Limitation

% -------------------------------------------------------------------------
%% ##### Eb/N0 Vector #####
% -------------------------------------------------------------------------
EbN0_file = open(f_getPath2file("EbN0",file_nbr));
EbN0 = EbN0_file.EbN0_vec;
EbN0_lin = 10.^(EbN0/10);
%
EbN0_steps = length(EbN0);

% -------------------------------------------------------------------------
%% ##### SNR and SNR (channel) Vector #####
% -------------------------------------------------------------------------
SNR_file = open(f_getPath2file("SNR",file_nbr));
SNR = SNR_file.SNR;
SNR_ch_file = open(f_getPath2file("SNR_ch",file_nbr));
SNR_ch = SNR_ch_file.SNR_ch;
%
N_ANT = size(SNR,1);
N_packetsPerEbN0 = size(SNR,2) / EbN0_steps;
%
if N_ANT > 0
   SNR_1 = SNR(1,:);
   SNR_1 = f_calc_meanSNRvec(SNR_1,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_1 = SNR_ch(1,:);
   SNR_ch_1 = f_calc_meanSNRvec(SNR_ch_1,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 1
   SNR_2 = SNR(2,:);
   SNR_2 = f_calc_meanSNRvec(SNR_2,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_2 = SNR_ch(2,:);
   SNR_ch_2 = f_calc_meanSNRvec(SNR_ch_2,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 2
   SNR_3 = SNR(3,:);
   SNR_3 = f_calc_meanSNRvec(SNR_3,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_3 = SNR_ch(3,:);
   SNR_ch_3 = f_calc_meanSNRvec(SNR_ch_3,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 3
   SNR_4 = SNR(4,:);
   SNR_4 = f_calc_meanSNRvec(SNR_4,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_4 = SNR_ch(4,:);
   SNR_ch_4 = f_calc_meanSNRvec(SNR_ch_4,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 4
   SNR_5 = SNR(5,:);
   SNR_5 = f_calc_meanSNRvec(SNR_5,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_5 = SNR_ch(5,:);
   SNR_ch_5 = f_calc_meanSNRvec(SNR_ch_5,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 5
   SNR_6 = SNR(6,:);
   SNR_6 = f_calc_meanSNRvec(SNR_6,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_6 = SNR_ch(6,:);
   SNR_ch_6 = f_calc_meanSNRvec(SNR_ch_6,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 6
   SNR_7 = SNR(7,:);
   SNR_7 = f_calc_meanSNRvec(SNR_7,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_7 = SNR_ch(7,:);
   SNR_ch_7 = f_calc_meanSNRvec(SNR_ch_7,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end
if N_ANT > 7
   SNR_8 = SNR(8,:);
   SNR_8 = f_calc_meanSNRvec(SNR_8,N_packetsPerEbN0,EbN0_steps,SNR_thres);
   SNR_ch_8 = SNR_ch(8,:);
   SNR_ch_8 = f_calc_meanSNRvec(SNR_ch_8,N_packetsPerEbN0,EbN0_steps,SNR_thres);
end

% -------------------------------------------------------------------------
%% ##### BER Vector #####
% -------------------------------------------------------------------------
BER_file = open(f_getPath2file("BER",file_nbr));
BER = BER_file.BER;
%
if N_ANT > 0
    BER_1 = BER(1,:);
    BER_1 = f_calc_meanBERvec(BER_1,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 1
    BER_2 = BER(2,:);
    BER_2 = f_calc_meanBERvec(BER_2,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 2
    BER_3 = BER(3,:);
    BER_3 = f_calc_meanBERvec(BER_3,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 3
    BER_4 = BER(4,:);
    BER_4 = f_calc_meanBERvec(BER_4,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 4
    BER_5 = BER(5,:);
    BER_5 = f_calc_meanBERvec(BER_5,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 5
    BER_6 = BER(6,:);
    BER_6 = f_calc_meanBERvec(BER_6,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 6
    BER_7 = BER(7,:);
    BER_7 = f_calc_meanBERvec(BER_7,N_packetsPerEbN0,EbN0_steps);
end
if N_ANT > 7
    BER_8 = BER(8,:);
    BER_8 = f_calc_meanBERvec(BER_8,N_packetsPerEbN0,EbN0_steps);
end

% -------------------------------------------------------------------------
%% ##### Ps and Ps (channel) Vector #####
% -------------------------------------------------------------------------
Ps_file = open(f_getPath2file("Ps",file_nbr));
Ps = Ps_file.Ps_vec;
Ps_ch_file = open(f_getPath2file("Ps_ch",file_nbr));
Ps_ch = Ps_ch_file.Ps_ch_vec;
%
N_packetsPerTX = EbN0_steps*N_packetsPerEbN0;
%
if N_ANT > 0
    Ps_1 = Ps(0*N_packetsPerTX+1 : 1*N_packetsPerTX);
    Ps_ch_1 = Ps_ch(0*N_packetsPerTX+1 : 1*N_packetsPerTX);
end
if N_ANT > 1
    Ps_2 = Ps(1*N_packetsPerTX+1 : 2*N_packetsPerTX);
    Ps_ch_2 = Ps_ch(1*N_packetsPerTX+1 : 2*N_packetsPerTX);
end
if N_ANT > 2
    Ps_3 = Ps(2*N_packetsPerTX+1 : 3*N_packetsPerTX);
    Ps_ch_3 = Ps_ch(2*N_packetsPerTX+1 : 3*N_packetsPerTX);
end
if N_ANT > 3
    Ps_4 = Ps(3*N_packetsPerTX+1 : 4*N_packetsPerTX);
    Ps_ch_4 = Ps_ch(3*N_packetsPerTX+1 : 4*N_packetsPerTX);
end
if N_ANT > 4
    Ps_5 = Ps(4*N_packetsPerTX+1 : 5*N_packetsPerTX);
    Ps_ch_5 = Ps_ch(4*N_packetsPerTX+1 : 5*N_packetsPerTX);
end
if N_ANT > 5
    Ps_6 = Ps(5*N_packetsPerTX+1 : 6*N_packetsPerTX);
    Ps_ch_6 = Ps_ch(5*N_packetsPerTX+1 : 6*N_packetsPerTX);
end
if N_ANT > 6
    Ps_7 = Ps(6*N_packetsPerTX+1 : 7*N_packetsPerTX);
    Ps_ch_7 = Ps_ch(6*N_packetsPerTX+1 : 7*N_packetsPerTX);
end
if N_ANT > 7
    Ps_8 = Ps(7*N_packetsPerTX+1 : 8*N_packetsPerTX);
    Ps_ch_8 = Ps_ch(7*N_packetsPerTX+1 : 8*N_packetsPerTX);
end

% -------------------------------------------------------------------------
%% ##### Pn and Pn (channel) Vector #####
% -------------------------------------------------------------------------
Pn_file = open(f_getPath2file("Pn",file_nbr));
Pn = Pn_file.Pn_vec;
Pn_ch_file = open(f_getPath2file("Pn_ch",file_nbr));
Pn_ch = Pn_ch_file.Pn_ch_vec;
%
if N_ANT > 0
    Pn_1 = Pn(0*N_packetsPerTX+1 : 1*N_packetsPerTX);
    Pn_ch_1 = Pn_ch(0*N_packetsPerTX+1 : 1*N_packetsPerTX);
end
if N_ANT > 1
    Pn_2 = Pn(1*N_packetsPerTX+1 : 2*N_packetsPerTX);
    Pn_ch_2 = Pn_ch(1*N_packetsPerTX+1 : 2*N_packetsPerTX);
end
if N_ANT > 2
    Pn_3 = Pn(2*N_packetsPerTX+1 : 3*N_packetsPerTX);
    Pn_ch_3 = Pn_ch(2*N_packetsPerTX+1 : 3*N_packetsPerTX);
end
if N_ANT > 3
    Pn_4 = Pn(3*N_packetsPerTX+1 : 4*N_packetsPerTX);
    Pn_ch_4 = Pn_ch(3*N_packetsPerTX+1 : 4*N_packetsPerTX);
end
if N_ANT > 4
    Pn_5 = Pn(4*N_packetsPerTX+1 : 5*N_packetsPerTX);
    Pn_ch_5 = Pn_ch(4*N_packetsPerTX+1 : 5*N_packetsPerTX);
end
if N_ANT > 5
    Pn_6 = Pn(5*N_packetsPerTX+1 : 6*N_packetsPerTX);
    Pn_ch_6 = Pn_ch(5*N_packetsPerTX+1 : 6*N_packetsPerTX);
end
if N_ANT > 6
    Pn_7 = Pn(6*N_packetsPerTX+1 : 7*N_packetsPerTX);
    Pn_ch_7 = Pn_ch(6*N_packetsPerTX+1 : 7*N_packetsPerTX);
end
if N_ANT > 7
    Pn_8 = Pn(7*N_packetsPerTX+1 : 8*N_packetsPerTX);
    Pn_ch_8 = Pn_ch(7*N_packetsPerTX+1 : 8*N_packetsPerTX);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% ##### 1 to 8 TX PLOTS ##### 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% ##### PLOT: SNR (RX) vs. EbN0 #####
% -------------------------------------------------------------------------
figure()
plot(EbN0,EbN0,usr_mrkr,'LineWidth',1.5*usr_LW,'MarkerSize',1.5*usr_MS);
hold on;
if N_ANT > 0
   plot(EbN0,SNR_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}"])
end
if N_ANT > 1
   plot(EbN0,SNR_2,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}"])
end
if N_ANT > 2
   plot(EbN0,SNR_3,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}"])
end
if N_ANT > 3
   plot(EbN0,SNR_4,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}"])
end
if N_ANT > 4
   plot(EbN0,SNR_5,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}"])
end
if N_ANT > 5
   plot(EbN0,SNR_6,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}"])
end
if N_ANT > 6
   plot(EbN0,SNR_7,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}";"SNR_{7TX}"])
end
if N_ANT > 7
   plot(EbN0,SNR_8,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}";"SNR_{7TX}";"SNR_{8TX}"])
end
hold off;
grid();
set(gca,'FontSize',usr_FS);
title("SNR (RX) vs. EbN0")
ylabel("SNR / EbN0 in dB")
xlabel("EbN0 in dB");
xlim([min(EbN0)-1; max(EbN0)+1])

% -------------------------------------------------------------------------
%% ##### PLOT: SNR (channel) vs. EbN0 #####
% -------------------------------------------------------------------------
figure()
plot(EbN0,EbN0,usr_mrkr,'LineWidth',1.5*usr_LW,'MarkerSize',1.5*usr_MS);
hold on;
if N_ANT > 0
   plot(EbN0,SNR_ch_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}"])
end
if N_ANT > 1
   plot(EbN0,SNR_ch_2,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}"])
end
if N_ANT > 2
   plot(EbN0,SNR_ch_3,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}"])
end
if N_ANT > 3
   plot(EbN0,SNR_ch_4,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}"])
end
if N_ANT > 4
   plot(EbN0,SNR_ch_5,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}"])
end
if N_ANT > 5
   plot(EbN0,SNR_ch_6,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}"])
end
if N_ANT > 6
   plot(EbN0,SNR_ch_7,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}";"SNR_{7TX}"])
end
if N_ANT > 7
   plot(EbN0,SNR_ch_8,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["EbN0";"SNR_{1TX}";"SNR_{2TX}";"SNR_{3TX}";"SNR_{4TX}";"SNR_{5TX}";"SNR_{6TX}";"SNR_{7TX}";"SNR_{8TX}"])
end
hold off;
grid();
set(gca,'FontSize',usr_FS);
title("SNR (Channel) vs. EbN0")
ylabel("SNR / EbN0 in dB")
xlabel("EbN0 in dB");
xlim([min(EbN0)-1; max(EbN0)+1])

% -------------------------------------------------------------------------
%% ##### PLOT: BER vs EbN0 #####
% -------------------------------------------------------------------------
figure()
if N_ANT > 0
   semilogy(EbN0,BER_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend("1 TX")
   hold on;
end
if N_ANT > 1
   semilogy(EbN0,BER_2,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX"])
end
if N_ANT > 2
   semilogy(EbN0,BER_3,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX"])
end
if N_ANT > 3
   semilogy(EbN0,BER_4,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX"])
end
if N_ANT > 4
   semilogy(EbN0,BER_5,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX"])
end
if N_ANT > 5
   semilogy(EbN0,BER_6,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX"])
end
if N_ANT > 6
   semilogy(EbN0,BER_7,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX"])
end
if N_ANT > 7
   semilogy(EbN0,BER_8,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX";"8 TX"])
end
hold off;
grid()
set(gca,'FontSize',usr_FS)
title("BER vs. EbN0")
ylabel("BER")
xlabel("EbN0 in dB");
xlim([min(EbN0)-1; max(EbN0)+1])

% -------------------------------------------------------------------------
%% ##### PLOT: BER vs SNR (expected) #####
% -------------------------------------------------------------------------
SNR_exp = EbN0 - 13;
% SNR_exp = EbN0 - 16;
%
figure()
if N_ANT > 0
   semilogy(SNR_exp,BER_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend("1 TX")
   hold on;
end
if N_ANT > 1
   semilogy(SNR_exp,BER_2,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX"])
end
if N_ANT > 2
   semilogy(SNR_exp,BER_3,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX"])
end
if N_ANT > 3
   semilogy(SNR_exp,BER_4,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX"])
end
if N_ANT > 4
   semilogy(SNR_exp,BER_5,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX"])
end
if N_ANT > 5
   semilogy(SNR_exp,BER_6,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX"])
end
if N_ANT > 6
   semilogy(SNR_exp,BER_7,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX"])
end
if N_ANT > 7
   semilogy(SNR_exp,BER_8,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
   legend(["1 TX";"2 TX";"3 TX";"4 TX";"5 TX";"6 TX";"7 TX";"8 TX"])
end
hold off;
grid()
set(gca,'FontSize',usr_FS)
title("BER vs. SNR")
ylabel("BER")
xlabel("SNR in dB");
xlim([min(SNR_exp)-1; max(SNR_exp)+1])


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% ##### 1 PLOTS ##### 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% ##### PLOT: SNR (RX) and SNR (channel) vs. EbN0 #####
% -------------------------------------------------------------------------
figure()
plot(EbN0,EbN0,usr_mrkr,'LineWidth',1.5*usr_LW,'MarkerSize',1.5*usr_MS); hold on;
plot(EbN0,SNR_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS);
plot(EbN0,SNR_ch_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS); hold off;
legend(["EbN0";"SNR_{RX}";"SNR_{Ch}"])
grid();
set(gca,'FontSize',usr_FS);
title("SNR (RX) and SNR (Ch) vs. EbN0")
ylabel("SNR / EbN0 in dB")
xlabel("EbN0 in dB");
xlim([min(EbN0)-1; max(EbN0)+1])

% -------------------------------------------------------------------------
%% ##### PLOT: Ps vs. Ps (channel) #####
% -------------------------------------------------------------------------
figure()
semilogy(Ps_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS); hold on;
semilogy(Ps_ch_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS); hold off;
legend(["P_{S,RX}";"P_{S,Ch}"])
grid();
set(gca,'FontSize',usr_FS);
title("P_{S,RX} and P_{S,Ch}")
ylabel("P_S")
xlabel("Transmission Index");

% -------------------------------------------------------------------------
%% ##### PLOT: Pn vs. Pn (channel) #####
% -------------------------------------------------------------------------
figure()
semilogy(Pn_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS); hold on;
semilogy(Pn_ch_1,usr_mrkr,'LineWidth',usr_LW,'MarkerSize',usr_MS); hold off;
legend(["P_{N,RX}";"P_{N,Ch}"])
grid();
set(gca,'FontSize',usr_FS);
title("P_{N,RX} and P_{N,Ch}")
ylabel("P_N")
xlabel("Transmission Index");