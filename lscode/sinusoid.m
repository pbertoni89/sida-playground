% Example: amplitude and phase of a sinusoid
F = 2;
omega = 2*pi*F;
T = [  2.188; 3.043;  4.207; 4.937;  5.675; 6.104;  6.260;  7.150; 8.600;  9.655 ];
Y = [ -1.112; 2.358; -1.807; 1.202; -0.814; 1.298; -2.520; -0.132; 1.421; -0.302 ];
Phi = [sin(omega*T), cos(omega*T)];
thetaLS = pinv(Phi)*Y;
Ahat    = sqrt(thetaLS(1)^2 + thetaLS(2)^2)
phihat  = atan2(thetaLS(2), thetaLS(1))
