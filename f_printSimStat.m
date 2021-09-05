function [] = f_printSimStat(currentRun,N_runs,EbN0_db,BER)

% Initialize status string
res = "Progress: Run " + num2str(currentRun);

% Add total number of runs if given
if nargin > 1
    res = res + " of " + num2str(N_runs);
end

% Add time stamp
res = res + ";\t Time: " + datestr(now,'HH:MM:SS');

% Add current EbN0 if given
if nargin > 2
    res = res + ";\t Current EbN0 = " + num2str(EbN0_db) + "dB";
end

% Add current BER if given
if nargin > 3
    res = res + ";\t Current BER = " + num2str(round(BER,11)*100) + " %%";
end

% Add line break
res = res + "\n";

% Print simulation status
fprintf(res);

end

