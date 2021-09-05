function [N_sc,N_g_l,N_g_r,N_pilot,N_data,bitsPerSymbol,N_bits,L_cp,...
            T_sym_OFDM,L_sym_OFDM,L_RX,pilot_idx,pilot_values,dc_idx] = f_getAnyParam_OFDM()
        
%F_GetAnyParam Returns all parameters to be used at any place.

%% General OFDM-Param

% Number of Subcarriers
N_sc = 64;

% Number of left Guard-Bands
N_g_l = 8;

% Number of right Guard-Bands
N_g_r = 7;

% Number of Pilot-Subcarrier
N_pilot = 8;

% Number of Payload-Subcarrier
N_data = N_sc - (N_g_l + N_g_r + N_pilot + 1); % DC not used

% Number of Payload-Bits per Subcarrier
bitsPerSymbol = 2;

% Number of Payload-Bits per OFDM-symbol
N_bits = N_data * bitsPerSymbol;


%% OFDM-Processing

% Length of Cyclic-Prefix (CP)
L_cp = 16;

% OFDM Symbol Duration
T_sym_OFDM = 1e-3;

% OFDM Symbol Length
L_sym_OFDM = N_sc + L_cp;

% OFDM symbol length with delay (padded zeros)
L_RX = 6e4;

% DC-index
dc_idx = N_sc/2+1;

% Pilot-Subcarrier Position
pilot_idx = [-24,-17,-10,-3,3,10,17,24] + dc_idx;

% Pilot-Subcarrier Value
pilot_values = 0 * ones(1, length(pilot_idx));

end

