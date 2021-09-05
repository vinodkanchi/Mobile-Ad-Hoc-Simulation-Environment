function u_tx5  = ...
    f_SIM_burst_tx5(ts1,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,N_TX,ts_scale)

%% Get necessary Param
[~,~,~,~,~,~,SPS,~,~,~,~,~,~,~,~,~,~] = f_getAnyParam_QAM();
[N_sc,~,~,~,~,~,~,~,T_sym_OFDM,~,~,~,~, ~] = f_getAnyParam_OFDM();

%% Scale for Eb = 1
c_d     = 9.927967304339253e+02;
c_ts    = 1.706666666666667e+03;
%
ts1 = sqrt(c_ts) * ts1;
% Avoid problems with synchronisation for low EbN0 (SNR)
if ts_scale < 1
    ts1 = sqrt(1/ts_scale) * ts1;  
end
%
data1 = sqrt(1/N_TX) * sqrt(c_d) * data1;
data2 = sqrt(1/N_TX) * sqrt(c_d) * data2;
data3 = sqrt(1/N_TX) * sqrt(c_d) * data3;
data4 = sqrt(1/N_TX) * sqrt(c_d) * data4;
data5 = sqrt(1/N_TX) * sqrt(c_d) * data5;
data6 = sqrt(1/N_TX) * sqrt(c_d) * data6;
data7 = sqrt(1/N_TX) * sqrt(c_d) * data7;
data8 = sqrt(1/N_TX) * sqrt(c_d) * data8;
data9 = sqrt(1/N_TX) * sqrt(c_d) * data9;
data10 = sqrt(1/N_TX) * sqrt(c_d) * data10;
data11 = sqrt(1/N_TX) * sqrt(c_d) * data11;
data12 = sqrt(1/N_TX) * sqrt(c_d) * data12;
data13 = sqrt(1/N_TX) * sqrt(c_d) * data13;
data14 = sqrt(1/N_TX) * sqrt(c_d) * data14;
data15 = sqrt(1/N_TX) * sqrt(c_d) * data15;
data16 = sqrt(1/N_TX) * sqrt(c_d) * data16;

%% Compose burst
b = [zeros(size(ts1));zeros(size(ts1));zeros(size(ts1));zeros(size(ts1));ts1;...
    zeros(size(ts1));zeros(size(ts1));zeros(size(ts1));...
        data1;data2;data3;data4;data5;data6;data7;data8;data9;...
        data10;data11;data12;data13;data14;data15;data16];

% Add imperfections (frequency offset; random phase added in channel)
dt = T_sym_OFDM/N_sc/SPS;
t = 0 : dt : dt*(length(b)-1);
FO = 0;
u_tx5 = b .* exp(1j*2*pi*FO*t.');

