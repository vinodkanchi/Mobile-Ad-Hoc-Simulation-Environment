function sym_decoded = f_SIM_OFDM_decode(sym_coded)

%% Get Param
[N_sc,N_g_l,N_g_r,~,~,~,~,L_cp,~,~,~,pilot_idx,~,dc_idx] = f_getAnyParam_OFDM();

%% Remove CP
sym_decoded = sym_coded(L_cp+1:end);

%% Calc FFT
sym_decoded = fftshift(fft(sym_decoded));

%% Remove Pilots
idx_g_l = 1 : 1 : N_g_l;
idx_g_r = N_sc-N_g_r+1 : 1 : N_sc;
idx_remove = [idx_g_l, pilot_idx, dc_idx, idx_g_r];
data_indeces = 1 : 1 : N_sc;
data_indeces(idx_remove) = [];
%
% sym_decoded = sym_decoded(data_indeces).';
sym_decoded = sym_decoded(data_indeces);

end

