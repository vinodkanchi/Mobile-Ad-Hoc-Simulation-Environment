function [t1,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,H5_ray] =... 
               f_SIM_OFDM_tx5(data1,data2,data3,data4,data5,data6,data7,data8,c_fading,c_fadeAllSC,c_singleCarrier)
           
%% Parameters
[N_sc,N_g_l,N_g_r,~,N_data,~,~,L_cp,~,~,~,pilot_idx,pilot_values,dc_idx]...
    = f_getAnyParam_OFDM();
[~,~,~,~,~,~,SPS,~,~,~,~,~,~,~,~,~,~] = f_getAnyParam_QAM();

%% Define Rayleigh Fading
H5_arr = complex(ones(N_data,1)); 
if c_fading
    if c_fadeAllSC
        H5_ray = sqrt(1/2)*(randn(N_data,1)+1j*(randn(N_data,1)));
    else
        H5_ray = sqrt(1/2)*(randn(1,1)+1j*(randn(1,1))) .* H5_arr;
    end
else
    H5_ray = H5_arr;
end
%% Define Single-Carrier Behaviour
data1_sc = data1(1);
data2_sc = data2(1);
data3_sc = data3(1);
data4_sc = data4(1);
data5_sc = data5(1);
data6_sc = data6(1);
data7_sc = data7(1);
data8_sc = data8(1);

%% ### Create Training Symbol ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, idx_g_r];

% Define all pilot index-values
pilot_indeces = 1 : 1 : N_sc;
pilot_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add payload information
y_tmp(pilot_indeces) = complex([...
    -1, -1,  1,  1, -1,  1, -1,  1,  1,  1, 1,  1,  1, -1, -1,  1, ...
    1, -1,  1, -1,  1,  1,  1,  1, 0, 1, -1, -1,  1,  1, -1,  1, -1,...
    1, -1, -1, -1, -1, -1, 1,  1, -1, -1,  1, -1,  1, -1,  1, 1].', 0);
      
%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Compose Training Symbol
t1 = [y_t(end/2+1:end); y_t; y_t; y_t(1:end/2)];
% t1 = [y_t(end/2+1:end); y_t; y_t; y_t; y_t; y_t; y_t; y_t(1:end/2)];



%% ### Create OFDM-Symbol 1 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = data5_sc;
else
    y_tmp(data_indeces) = data5(1:length(data_indeces));
end

% Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y1 = [y_t(end-(SPS*L_cp-1):end); y_t];


%% ### Create OFDM-Symbol 2 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = data6_sc;
else
    y_tmp(data_indeces) = data6(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y2 = [y_t(end-(SPS*L_cp-1):end); y_t];
%% ### Create OFDM-Symbol 3 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = data7_sc;
else
    y_tmp(data_indeces) = data7(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y3 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 4 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = data8_sc;
else
    y_tmp(data_indeces) = data8(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;
%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y4 = [y_t(end-(SPS*L_cp-1):end); y_t];


%% ### Create OFDM-Symbol 5 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = data1_sc;
else
    y_tmp(data_indeces) = data1(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;
%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y5 = [y_t(end-(SPS*L_cp-1):end); y_t];
%% ### Create OFDM-Symbol 6 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = -data2_sc;
else
    y_tmp(data_indeces) = -data2(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y6 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 7 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = -data3_sc;
else
    y_tmp(data_indeces) = -data3(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y7 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 8 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = -data4_sc;
else
    y_tmp(data_indeces) = -data4(1:length(data_indeces));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y8 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 9 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

if c_singleCarrier
    y_tmp(dc_idx+1) = conj(data5_sc);
else
    y_tmp(data_indeces) = conj(data5(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y9 = [y_t(end-(SPS*L_cp-1):end); y_t];
%% ### Create OFDM-Symbol 10 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = conj(data6_sc);
else
    y_tmp(data_indeces) = conj(data6(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y10 = [y_t(end-(SPS*L_cp-1):end); y_t];
%% ### Create OFDM-Symbol 11 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = conj(data7_sc);
else
    y_tmp(data_indeces) = conj(data7(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y11 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 12 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = conj(data8_sc);
else
    y_tmp(data_indeces) = conj(data8(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y12 = [y_t(end-(SPS*L_cp-1):end); y_t];


%% ### Create OFDM-Symbol 13 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = conj(data1_sc);
else
    y_tmp(data_indeces) = conj(data1(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y13 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 14 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = -conj(data2_sc);
else
    y_tmp(data_indeces) = -conj(data2(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y14 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 15 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = -conj(data3_sc);
else
    y_tmp(data_indeces) = -conj(data3(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y15 = [y_t(end-(SPS*L_cp-1):end); y_t];

%% ### Create OFDM-Symbol 16 ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Allocate memory
y_tmp = complex(zeros(N_sc, 1));

% Add pilot information
y_tmp(pilot_idx) = pilot_values;

% Add payload information (--> According to Code Matrix)
if c_singleCarrier
    y_tmp(dc_idx+1) = -conj(data4_sc);
else
    y_tmp(data_indeces) = -conj(data4(1:length(data_indeces)));
end
 % Perform Fading in Frequency Domain
y_tmp(data_indeces) = y_tmp(data_indeces) .* H5_ray;

%% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

%% Perform IFFT
% Use ifftshift to avoid problems when calling the MATLAB function ifft
% ifftshift shifts the DC components to the middle of the spectrum
y_t = ifft(ifftshift(y_tmp));

%% Add CP
y16 = [y_t(end-(SPS*L_cp-1):end); y_t];
