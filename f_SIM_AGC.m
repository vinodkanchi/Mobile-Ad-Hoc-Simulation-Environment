function y = f_SIM_AGC(u,c_AGC)

%% Get necessary param
% [~,~,viewPlots,RX_const_type,~,~] = f_getAnyParam_Simulation();

%% Perform AGC
if c_AGC == 0
    y = u;
else
    % Turn all symbols to the same angle
    y_sameAngle = u .^ 16;
    
    % Calc mean of absolute values
    y_mean = mean(abs(y_sameAngle));
    
    % Calc root again as factor for correction
    y_correct = y_mean ^ (1/16);
    
    % Normalize
    y = u / y_correct;
end


%% Plot Constellation
% % if viewPlots
% %     figure(5)
% %     subplot(2,2,1)
% %     if RX_const_type == "scatter"
% %         scatter(real(y), imag(y), 'filled')
% %     elseif RX_const_type == "plot"
% %         plot(y)
% %     end
% %     grid()
% %     title("Constellation: Symbol 1");
% %     xlabel("\Re(sym1)"); ylabel("\Im(sym1)")
% %     axis([-1 1 -1 1] * 1.3);
% % end
