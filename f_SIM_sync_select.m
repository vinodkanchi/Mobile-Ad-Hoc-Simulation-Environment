function [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,...
    H1e,H2e,H3e,H4e,H5e,H6e,H7e,H8e,SNR] = ...
    f_SIM_sync_select(y11,y12,y13,y14,y15,y16,y17,y18,y19,y110,y111,y112,y113,y114,y115,y116,...
                H11e,H12e,H13e,H14e,H15e,H16e,H17e,H18e,...
                y21,y22,y23,y24,y25,y26,y27,y28,y29,y210,y211,y212,y213,y214,y215,y216,...
                H21e,H22e,H23e,H24e,H25e,H26e,H27e,H28e,...
                y31,y32,y33,y34,y35,y36,y37,y38,y39,y310,y311,y312,y313,y314,y315,y316,...
                H31e,H32e,H33e,H34e,H35e,H36e,H37e,H38e,...
                y41,y42,y43,y44,y45,y46,y47,y48,y49,y410,y411,y412,y413,y414,y415,y416,...
                H41e,H42e,H43e,H44e,H45e,H46e,H47e,H48e,...
                y51,y52,y53,y54,y55,y56,y57,y58,y59,y510,y511,y512,y513,y514,y515,y516,...
                H51e,H52e,H53e,H54e,H55e,H56e,H57e,H58e,...
                y61,y62,y63,y64,y65,y66,y67,y68,y69,y610,y611,y612,y613,y614,y615,y616,...
                H61e,H62e,H63e,H64e,H65e,H66e,H67e,H68e,...
                y71,y72,y73,y74,y75,y76,y77,y78,y79,y710,y711,y712,y713,y714,y715,y716,...
                H71e,H72e,H73e,H74e,H75e,H76e,H77e,H78e,...
                y81,y82,y83,y84,y85,y86,y87,y88,y89,y810,y811,y812,y813,y814,y815,y816,...
                H81e,H82e,H83e,H84e,H85e,H86e,H87e,H88e,...
                Ps1,Pn1,Ps2,Pn2,Ps3,Pn3,Ps4,Pn4,Ps5,Pn5,Ps6,Pn6,Ps7,Pn7,Ps8,Pn8)

H1 = sum(abs(H11e.^2));
H2 = sum(abs(H22e.^2));
H3 = sum(abs(H33e.^2));
H4 = sum(abs(H44e.^2));
H5 = sum(abs(H55e.^2));
H6 = sum(abs(H66e.^2));
H7 = sum(abs(H77e.^2));
H8 = sum(abs(H88e.^2));
% H1 = sum(abs((H11e.*(0.9963 + 1.0021j)).^2));
% H2 = sum(abs((H22e.*(0.4748 - 0.8538j)).^2));
% H3 = sum(abs((H33e.*(0.5140 - 0.2146j)).^2));
% H4 = sum(abs((H44e.*(-2.0819 + 1.0171j)).^2));
            
if (H1 > H2) && (H1 > H3) && (H1 > H4) && (H1 > H5) && (H1 > H6) && (H1 > H7) && (H1 > H8)
    d1 = y11;
    d2 = y12;
    d3 = y13;
    d4 = y14;
    d5 = y15;
    d6 = y16;
    d7 = y17;
    d8 = y18;
    d9 = y19;
    d10 = y110;
    d11 = y111;
    d12 = y112;
    d13 = y113;
    d14 = y114;
    d15 = y115;
    d16 = y116;
    H1e = H11e;
    H2e = H12e;
    H3e = H13e;
    H4e = H14e;
    H5e = H15e;
    H6e = H16e;
    H7e = H17e;
    H8e = H18e;
    SNR = 10*log10(Ps1/Pn1);
elseif (H2 > H1) && (H2 > H3) && (H2 > H4) && (H2 > H5) && (H2 > H6) && (H2 > H7) && (H2 > H8)
    d1 = y21;
    d2 = y22;
    d3 = y23;
    d4 = y24;
    d5 = y25;
    d6 = y26;
    d7 = y27;
    d8 = y28;
    d9 = y29;
    d10 = y210;
    d11 = y211;
    d12 = y212;
    d13 = y213;
    d14 = y214;
    d15 = y215;
    d16 = y216;
    H1e = H21e;
    H2e = H22e;
    H3e = H23e;
    H4e = H24e;
    H5e = H25e;
    H6e = H26e;
    H7e = H27e;
    H8e = H28e;
    SNR = 10*log10(Ps2/Pn2);
elseif (H3 > H1) && (H3 > H2) && (H3 > H4) && (H3 > H5) && (H3 > H6) && (H3 > H7) && (H3 > H8)
    d1 = y31;
    d2 = y32;
    d3 = y33;
    d4 = y34;
    d5 = y35;
    d6 = y36;
    d7 = y37;
    d8 = y38;
    d9 = y39;
    d10 = y310;
    d11 = y311;
    d12 = y312;
    d13 = y313;
    d14 = y314;
    d15 = y315;
    d16 = y316;
    H1e = H31e;
    H2e = H32e;
    H3e = H33e;
    H4e = H34e;
    H5e = H35e;
    H6e = H36e;
    H7e = H37e;
    H8e = H38e;
    SNR = 10*log10(Ps3/Pn3);
elseif (H4 > H1) && (H4 > H2) && (H4 > H3) && (H4 > H5) && (H4 > H6) && (H4 > H7) && (H4 > H8)
    d1 = y41;
    d2 = y42;
    d3 = y43;
    d4 = y44;
    d5 = y45;
    d6 = y46;
    d7 = y47;
    d8 = y48;
    d9 = y49;
    d10 = y410;
    d11 = y411;
    d12 = y412;
    d13 = y413;
    d14 = y414;
    d15 = y415;
    d16 = y416;
    H1e = H41e;
    H2e = H42e;
    H3e = H43e;
    H4e = H44e;
    H5e = H45e;
    H6e = H46e;
    H7e = H47e;
    H8e = H48e;
    SNR = 10*log10(Ps4/Pn4);
elseif (H5 > H1) && (H5 > H2) && (H5 > H3) && (H5 > H4) && (H5 > H6) && (H5 > H7) && (H5 > H8)
    d1 = y51;
    d2 = y52;
    d3 = y53;
    d4 = y54;
    d5 = y55;
    d6 = y56;
    d7 = y57;
    d8 = y58;
    d9 = y59;
    d10 = y510;
    d11 = y511;
    d12 = y512;
    d13 = y513;
    d14 = y514;
    d15 = y515;
    d16 = y516;
    H1e = H51e;
    H2e = H52e;
    H3e = H53e;
    H4e = H54e;
    H5e = H55e;
    H6e = H56e;
    H7e = H57e;
    H8e = H58e;
    SNR = 10*log10(Ps5/Pn5);
elseif (H6 > H1) && (H6 > H2) && (H6 > H3) && (H6 > H4) && (H6 > H5) && (H6 > H7) && (H6 > H8)
    d1 = y61;
    d2 = y62;
    d3 = y63;
    d4 = y64;
    d5 = y65;
    d6 = y66;
    d7 = y67;
    d8 = y68;
    d9 = y69;
    d10 = y610;
    d11 = y611;
    d12 = y612;
    d13 = y613;
    d14 = y614;
    d15 = y615;
    d16 = y616;
    H1e = H61e;
    H2e = H62e;
    H3e = H63e;
    H4e = H64e;
    H5e = H65e;
    H6e = H66e;
    H7e = H67e;
    H8e = H68e;
    SNR = 10*log10(Ps6/Pn6);
elseif (H7 > H1) && (H7 > H2) && (H7 > H3) && (H7 > H4) && (H7 > H5) && (H7 > H6) && (H7 > H8)
    d1 = y71;
    d2 = y72;
    d3 = y73;
    d4 = y74;
    d5 = y75;
    d6 = y76;
    d7 = y77;
    d8 = y78;
    d9 = y79;
    d10 = y710;
    d11 = y711;
    d12 = y712;
    d13 = y713;
    d14 = y714;
    d15 = y715;
    d16 = y716;
    H1e = H71e;
    H2e = H72e;
    H3e = H73e;
    H4e = H74e;
    H5e = H75e;
    H6e = H76e;
    H7e = H77e;
    H8e = H78e;
    SNR = 10*log10(Ps7/Pn7);
else 
    d1 = y81;
    d2 = y82;
    d3 = y83;
    d4 = y84;
    d5 = y85;
    d6 = y86;
    d7 = y87;
    d8 = y88;
    d9 = y89;
    d10 = y810;
    d11 = y811;
    d12 = y812;
    d13 = y813;
    d14 = y814;
    d15 = y815;
    d16 = y816;
    H1e = H81e;
    H2e = H82e;
    H3e = H83e;
    H4e = H84e;
    H5e = H85e;
    H6e = H86e;
    H7e = H87e;
    H8e = H88e;
    SNR = 10*log10(Ps8/Pn8);
    
end