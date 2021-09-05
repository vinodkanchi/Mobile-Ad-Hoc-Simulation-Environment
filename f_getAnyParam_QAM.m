function [fc,fi,Tf,Tsym,Bpf,Spf,SPS,L0,Ltx,MCR,ID,Ts,Tg,fsample,M,Lm,Lrx] ...
            = f_getAnyParam_QAM()
%F_GetAnyParam Return all parameters to be used at any place.

%%

% IF Frequency
fi = 20000;

% Carrier Frequency
fc = 80e6;
% fc = 100e6;
% fc = 1e9;

% Wertigkeit der QAM - Do not change, adaptation of code required
M = 4;

% Frame period [s]
Tf = 1; 

% Symbol rate [Hz]
Rs = 5000;

% Bits per Frame
Bpf = 255*log2(M);

% Samples per Symbol
SPS = 20;

% Length of m-sequence
Lm  = 31; 

% Length of receive sequence
Lrx = 10000; 

% Number of Symbols: Bits per Frame (with header/tail); do not change
Spf  = (Lm + 2) * 3 + Bpf/log2(M);

% Number of zeros prepended/appended to transmit sequence
L0 = (SPS * Tf* Rs - Spf * SPS);  

% Length of transmit sequence with zeros
Ltx  = L0 + Spf*SPS;

% Master Clock rate
MCR = 5e6;

% interpolation/decimation factor
ID = MCR*Tf/Ltx;

% Sample time
Ts = ID/MCR;

% Symbol duration
Tsym = 1/Rs;

% Duration of one transmit bit (adjusted in source)
Tg = Tf/Bpf;

% Sample frequency on SDR input/output
fsample = SPS / Tsym;

end

