function meanBERvec = f_calc_meanBERvec(BERvec_ext,SNR_range,SNR_steps)

meanBERvec = zeros(1, SNR_steps);

for i = 1 : 1 : SNR_steps
    meanBERvec(i) = mean(BERvec_ext((i-1)*SNR_range+1 : i*SNR_range));
end

end

