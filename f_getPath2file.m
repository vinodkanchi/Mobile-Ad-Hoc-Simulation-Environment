function path2file = f_getPath2file(fileType,fileNbr)

% -------------------------------------------------------------------------
%% BER Files
% -------------------------------------------------------------------------
BER_files = [...
    "./SIM_results/resfile_BER.mat";...
    ""  ];

% -------------------------------------------------------------------------
%% EbN0 Files
% -------------------------------------------------------------------------
EbN0_files = [...
    "./SIM_results/resfile_EbN0.mat";...
    ""  ];

% -------------------------------------------------------------------------
%% SNR Files
% -------------------------------------------------------------------------
SNR_files = [...
    "./SIM_results/resfile_SNR.mat";...
    ""  ];

SNR_ch_files = [...
    "./SIM_results/resfile_SNR_ch.mat";...
    ""  ];

% -------------------------------------------------------------------------
%% Ps Files
% -------------------------------------------------------------------------
Ps_files = [...
    "./SIM_results/resfile_Ps.mat";...
    ""  ];

Ps_ch_files = [...
    "./SIM_results/resfile_Ps_ch.mat";...
    ""  ];

% -------------------------------------------------------------------------
%% Pn Files
% -------------------------------------------------------------------------
Pn_files = [...
    "./SIM_results/resfile_Pn.mat";...
    ""  ];

Pn_ch_files = [...
    "./SIM_results/resfile_Pn_ch.mat";...
    ""  ];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Return Path to File
if fileType == "BER"
    path2file = BER_files(fileNbr);
elseif fileType == "EbN0"
    path2file = EbN0_files(fileNbr);
elseif fileType == "SNR"
    path2file = SNR_files(fileNbr);
elseif fileType == "SNR_ch"
    path2file = SNR_ch_files(fileNbr);
elseif fileType == "Ps"
    path2file = Ps_files(fileNbr);    
elseif fileType == "Ps_ch"
    path2file = Ps_ch_files(fileNbr);
elseif fileType == "Pn"
    path2file = Pn_files(fileNbr);
elseif fileType == "Pn_ch"
    path2file = Pn_ch_files(fileNbr);
else
    path2file = "";
end


end

