function [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12s,y13,y14,y15,y16,Ps7,Pn7] = ...
         f_SIM_sync_tx7(rx,sync_thr,sync_area,c_filt,c_sync,rnd_pos)

%% Get Param
coder.extrinsic('legend');
[N_sc,N_g_l,N_g_r,~,N_data,~,~,~,T_sym_OFDM,~,~,pilot_idx,~,dc_idx] = f_getAnyParam_OFDM();
[~,~,~,~,~,~,SPS,~,~,~,~,~,~,~,~,~,~] = f_getAnyParam_QAM();
% [~,~,viewPlots,~,~,~] = f_getAnyParam_Simulation();

sym_length = 2*N_data;

SNR_start = 53;


%% ### Create Training Symbol (Preamble) ###
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
      
% Downsampled preamble
y_t_ds = ifft(ifftshift(y_tmp));
preamble_ds = [y_t_ds(end/2+1:end); y_t_ds; y_t_ds; y_t_ds(1:end/2)];
ts_length = length(preamble_ds);
      
% Perform Upsampling in Frequency Domain (--> Sinc-Interpolation)
y_tmp = [zeros(1/2*(SPS-1)*length(y_tmp), 1);...
            y_tmp;...
            zeros(1/2*(SPS-1)*length(y_tmp), 1)];

% Perform IFFT
y_t = ifft(ifftshift(y_tmp));

% Compose upsampled preamble (--> used for correlation)
preamble = [y_t(end/2+1:end); y_t; y_t; y_t(1:end/2)];


%% ### Filter input signal ###
N_filt = 64;
% fs = 1; fpass = 0.055; fstop = 0.100
if (c_filt == 1)
    b_fir_ofdm = [-0.00067908;0.0010149;0.0013273;0.0016928;0.0017751;0.001336;0.00027722;-0.0013036;-0.0030899;-0.0045845;-0.0052167;-0.0045037;-0.0022289;0.0014142;0.0057648;0.0097624;0.012153;0.011788;0.0079892;0.00084347;-0.0086197;-0.018464;-0.026111;-0.028806;-0.024208;-0.010969;0.010829;0.039492;0.071855;0.10375;0.13068;0.14868;0.155;0.14868;0.13068;0.10375;0.071855;0.039492;0.010829;-0.010969;-0.024208;-0.028806;-0.026111;-0.018464;-0.0086197;0.00084347;0.0079892;0.011788;0.012153;0.0097624;0.0057648;0.0014142;-0.0022289;-0.0045037;-0.0052167;-0.0045845;-0.0030899;-0.0013036;0.00027722;0.001336;0.0017751;0.0016928;0.0013273;0.0010149;-0.00067908];
else
    b_fir_ofdm = [zeros(N_filt-1,1);1];
end
y_tmp = conv(rx, b_fir_ofdm);
y = y_tmp(N_filt/2+1 : N_filt/2+length(rx)); % y is the received filtered Signal.

% % Plot spectrum before downsampling --> plotted already in sync_tx1
% figure(4)
% subplot(2, 2, 4)
% [curr_spectrum, freq_norm] = ...
%     periodogram(y, hann(length(y), 'periodic'), length(y)*2, 1);
% curr_psd = 10*log10(fftshift(curr_spectrum));
% freq_norm = freq_norm - 0.5;
% plot(freq_norm, curr_psd); grid();
% xlim([-0.1, 0.1]); ylim([-35, +35]);
% title("Spectrum RX4 before Downsampling")


%% ### Sync ###

%% STEP 1
% Finding the convolution matrix XC.
% Correlating with second half of preamble
XC = (conv(y.',fliplr(preamble')));


%% STEP 2
% STEP 2a: Find the peak in the correlation result
[maxim,Index] = max(abs(XC(:)));

% STEP 2b: Downsampling
start_index = mod(Index, SPS) + 1; % verified w. simulations
y_ds = y(start_index:SPS:end);

% STEP 2c: Find the peak in the downsampled correlation result
XC_ds = (conv(y_ds.',fliplr(preamble_ds'))); % Correlating with second 
% We assume TDMA and we know where to search +/- 20 samples
[maxim,Index] = ...
    max(abs(XC_ds(7*ts_length-sync_area:7*ts_length+sync_area)));
if (c_sync == 0)
    Index = round(rnd_pos/SPS) + 7*ts_length + 1;
else
    Index = Index-sync_area + 7*ts_length - 1;
end

% Visualization
% % if viewPlots
% %     figure(2)
% %     subplot(4, 4, 5)
% %     plot(abs(XC_ds)); hold on
% %     plot(Index, abs(XC_ds(Index)), 'ro'); hold off; grid on
% %     title("Sync Tx5");
% % end


%% STEP 3
% Extract the required amount of samples from XC.
% Feed them as rows by shifting, into Matrix A.
% Each row of A is a window of length:'samples'.
samples = 17; % Total No.of samples
A = complex(zeros(samples,samples));
for k = 1:samples
   A(k,:) = XC_ds(Index+k-1-(samples-1):Index+k-1);
end


%% STEP 4
% Find the window with maximum power.
pow = zeros(1,samples); % Contains power values of each window
power = 0;
for i = 1:samples
   for j = 1:samples
       power = power+abs(A(i,j))^2;
   end
   pow(i) = power;
   power = 0;
end
[~,WindowNum] = max(pow(:));    % Finding the window with maximum power              
window = abs(A(WindowNum,:));   % Extract the window from A matrix.


%% STEP 5
% Find the first peak passing the threshold.
[maximum,indexofmax] = max(abs(window(:)));
check = 0;
Reqpeak = 0;        
for i = 1:samples
   if window(i)>(sync_thr*maximum) % Threshold found by trial and error
       Reqpeak = i;
       check = 1;
   end
   if check == 1
       break;
   end
end


%% STEP 6
% Output the detected beginning of burst
pv = zeros(1,1);
pv = Index - indexofmax + Reqpeak - 7*length(preamble)/SPS;
if pv < 0
    pv = 0;
end


%% ### Perform Channel and Frequency Offset Estimation ###
% Define all non-data index-values
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];

% Define all data index-values
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];

% Channel Estimation
rx_preamble_ds1 = y_ds(pv+0*ts_length+33 : pv+0*ts_length+33+63);
rx_preamble_ds2 = y_ds(pv+1*ts_length+33 : pv+1*ts_length+33+63);
rx_preamble_ds3 = y_ds(pv+2*ts_length+33 : pv+2*ts_length+33+63);
rx_preamble_ds4 = y_ds(pv+3*ts_length+33 : pv+3*ts_length+33+63);
rx_preamble_ds5 = y_ds(pv+4*ts_length+33 : pv+4*ts_length+33+63);
rx_preamble_ds6 = y_ds(pv+5*ts_length+33 : pv+5*ts_length+33+63);
rx_preamble_ds7 = y_ds(pv+6*ts_length+33 : pv+6*ts_length+33+63);
rx_preamble_ds8 = y_ds(pv+7*ts_length+33 : pv+7*ts_length+33+63);

RX_preamble_ds1 = fftshift(fft(rx_preamble_ds1));
RX_preamble_ds2 = fftshift(fft(rx_preamble_ds2));
RX_preamble_ds3 = fftshift(fft(rx_preamble_ds3));
RX_preamble_ds4 = fftshift(fft(rx_preamble_ds4));
RX_preamble_ds5 = fftshift(fft(rx_preamble_ds5));
RX_preamble_ds6 = fftshift(fft(rx_preamble_ds6));
RX_preamble_ds7 = fftshift(fft(rx_preamble_ds7));
RX_preamble_ds8 = fftshift(fft(rx_preamble_ds8));

Preamble_ds     = fftshift(fft(y_t_ds));

H71e_tmp = RX_preamble_ds1(data_indeces) ./ Preamble_ds(data_indeces);
H72e_tmp = RX_preamble_ds2(data_indeces) ./ Preamble_ds(data_indeces);
H73e_tmp = RX_preamble_ds3(data_indeces) ./ Preamble_ds(data_indeces);
H74e_tmp = RX_preamble_ds4(data_indeces) ./ Preamble_ds(data_indeces);
H75e_tmp = RX_preamble_ds5(data_indeces) ./ Preamble_ds(data_indeces);
H76e_tmp = RX_preamble_ds6(data_indeces) ./ Preamble_ds(data_indeces);
H77e_tmp = RX_preamble_ds7(data_indeces) ./ Preamble_ds(data_indeces);
H78e_tmp = RX_preamble_ds8(data_indeces) ./ Preamble_ds(data_indeces);

% Necessary due to Simulink's fixed data propagation
H71e = complex(zeros(N_data,1)); 
H72e = complex(zeros(N_data,1));
H73e = complex(zeros(N_data,1));
H74e = complex(zeros(N_data,1));
H75e = complex(zeros(N_data,1));
H76e = complex(zeros(N_data,1));
H77e = complex(zeros(N_data,1));
H78e = complex(zeros(N_data,1));
H71e = H71e + H71e_tmp;
H72e = H72e + H72e_tmp;
H73e = H73e + H73e_tmp;
H74e = H74e + H74e_tmp;
H75e = H75e + H75e_tmp;
H76e = H76e + H76e_tmp;
H77e = H77e + H77e_tmp;
H78e = H78e + H78e_tmp;

H7e = H71e(1:N_data) + H72e(1:N_data) + H73e(1:N_data) + H74e(1:N_data) + H75e(1:N_data) + H76e(1:N_data) + H77e(1:N_data) + H78e(1:N_data);
H7e = H7e / 8;

% Estimate Frequency Offset
FO7_est = angle(H71e'*H78e)/(2*pi*T_sym_OFDM);

% Visualization
% % if viewPlots
% %     figure(3)
% %     subplot(2,2,4); hold off
% %     plot(abs(H41e)); hold on
% %     plot(abs(H42e));
% %     plot(abs(H43e));
% %     plot(abs(H44e));
% %     plot(abs(H4e),'g')
% %     legend('H41e','H42e','H43e','H44e','H4e'); title('H4e(f)')
% %     xlim([0, 40]);
% % %     axis([0,40, 0,2]);
% %     grid();
% % end


%% Get 16 OFDM Symbols after TS4
y12 = y_ds(pv+1+8*192 : pv+8*192+16*sym_length);
y1 = y12(0*sym_length+1 : 1*sym_length);
y2 = y12(1*sym_length+1 : 2*sym_length);
y3 = y12(2*sym_length+1 : 3*sym_length);
y4 = y12(3*sym_length+1 : 4*sym_length);
y5 = y12(4*sym_length+1 : 5*sym_length);
y6 = y12(5*sym_length+1 : 6*sym_length);
y7 = y12(6*sym_length+1 : 7*sym_length);
y8 = y12(7*sym_length+1 : 8*sym_length);
y9 = y12(8*sym_length+1 : 9*sym_length);
y10 = y12(9*sym_length+1 : 10*sym_length);
y11 = y12(10*sym_length+1 : 11*sym_length);
y12s = y12(11*sym_length+1 : 12*sym_length);
y13 = y12(12*sym_length+1 : 13*sym_length);
y14 = y12(13*sym_length+1 : 14*sym_length);
y15 = y12(14*sym_length+1 : 15*sym_length);
y16 = y12(15*sym_length+1 : 16*sym_length);


%% Estimate SNR
% Ps7 = mean(abs(y12(3:end)).^2);
Ps7 = mean(abs(y12(SNR_start:end)).^2);
Pn7 = mean(abs(y_ds(pv+8*192+16*sym_length+1:end)).^2);

Ps7 = Ps7 - Pn7;
if Ps7 <= 0
    Ps7 = eps; % avoid negative power
end


end