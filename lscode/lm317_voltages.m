%---------------------------------------------------------
% Interpolation of the dependence between temperature and reference
% voltage of the LM317 integrated circuit
% Datasheet: http://www.ti.com.cn/cn/lit/ds/symlink/lm117.pdf
% Code: courtesy of Patrizio Bertoni, March 8 2014
%---------------------------------------------------------

format long e;    % Show enough significant digits


% Explanatory variable: temperatures from -50 to 125 degrees Celsius
temp = (-50:25:125)';

% Variable to be explained: reference voltages read empirically 
% (with pen and ruler) from page 8 of the datasheet, figure 1
voltages = [1.245 1.248 1.249 1.250 1.249 1.246 1.243 1.237]';

% Interpolation with a 2nd-degree polynomial, y = a_0 + a_1 t + a_2 t^2
onesv = ones(length(temp),1);
temp2 = temp.^2;
Phi2  = [onesv, temp, temp2];
thetaLS2 = pinv(Phi2)*voltages
% Notice that the third coefficient, although small, is negative:
% the parabola is concave

% Interpolation with a 3rd-degree polynomial, y = a_0 + a_1 t + a_2 t^2 + a_3 t^3
temp3 = temp.^3;
Phi3  = [Phi2, temp3];
thetaLS3 = pinv(Phi3)*voltages

% Example: "prediction" (extrapolation) of the ref. voltage at 150 deg. Celsius
new_temp = 150;
regressor2    = [1; new_temp; new_temp^2];
volt_predict2 = regressor2' * thetaLS2     % \varphi_{N+1}^\top \theta
regressor3    = [regressor2; new_temp^3];
volt_predict3 = regressor3' * thetaLS3     % \varphi_{N+1}^\top \theta


%---------------------------------------------------------
% Plot measures, 2nd-order model, and 3rd-order model
%---------------------------------------------------------
t = (-50:5:150)';
figure(1);
clf;
hold on;
grid on;

plot(temp, voltages, 'k*', 'LineStyle','none')
plot_volt2 = [ones(length(t),1), t, t.^2]*thetaLS2;
plot(t, plot_volt2, "r");
plot_volt3 = [ones(length(t),1), t, t.^2, t.^3]*thetaLS3;
plot(t, plot_volt3, "b");

axis([-50 150 1.225 1.255]);
title('LM317 temperature stability');
xlabel('Temperature');
ylabel('Ref. voltage');
legend('measures', '2nd-order model', '3rd-order model');
pause
