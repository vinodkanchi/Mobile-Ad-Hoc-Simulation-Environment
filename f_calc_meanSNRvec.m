function meanSNRvec = f_calc_meanSNRvec(SNRvec_ext,SNR_range,SNR_steps,SNR_threshold)

meanSNRvec = zeros(1, SNR_steps); 

for i = 1 : 1 : SNR_steps
    % Get SNR values of corresponding range
    SNR_tmp_vec =  SNRvec_ext((i-1)*SNR_range+1 : i*SNR_range);
    % Discard all SNR values below a user defined threshold
    SNR_aboveThres = SNR_tmp_vec > SNR_threshold;
    SNR_tmp_vec = SNR_tmp_vec .* SNR_aboveThres;
    % Calc mean value
    SNR_tmp_vec = 10.^(SNR_tmp_vec./10);
    SNR_tmp_vec = SNR_tmp_vec .* (SNR_tmp_vec ~= 1);
    meanSNRvec(i) = 10*log10(sum(SNR_tmp_vec)/sum(SNR_aboveThres));
end

end

