function [res] = f_convertElapsedTime(t_elapsed)
% Converts measured elapsed time in seconds into a readable string for
% console information.

if (t_elapsed > 60*60)
    t_h     = floor(t_elapsed / (60*60));
    t_min   = floor((t_elapsed-t_h*60*60) / 60);
    t_sec   = floor(rem(t_elapsed,60));
    res     = num2str(t_h) + " h " + num2str(t_min) + " min " + ...
                num2str(t_sec) + " sec";

elseif (t_elapsed > 60)
    t_min   = floor(t_elapsed / 60);
    t_sec   = floor(rem(t_elapsed,60));
    res     = num2str(t_min) + " min " + num2str(t_sec) + " sec";

else
    t_sec   = t_elapsed;
    res     = num2str(t_sec) + " sec";
end

end

