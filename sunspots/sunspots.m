% Octave/Matlab code:
% Autoregressive fit of Wolf's sunspot numbers
% Comparison of the power spectral densities estimates obtained from DFT and AR fit 
clc, clear all
close all
%-------------------------------------------------------------------
% Load data
sunspot_data = load('sunspots.txt');
N = size(sunspot_data,1);
yrs = sunspot_data(:, 1);     % Years
ss  = sunspot_data(:, 2);     % Wolf's sunspot numbers
ssd = ss - mean(ss);          % Detrend the time series

%-------------------------------------------------------------------
% Plot sunspot numbers
figure('name', 'Sunspots'), title('Sunspots'), clf(), plot(yrs, ss, 'r')
%hold on, plot(yrs, sin(0.56*yrs)+mean(ss), 'b')
axis([1749, 1927, 0, 160]), xlabel('year'), ylabel('Wolf''s number')

%-------------------------------------------------------------------
% Estimate the power spectrum by PEM (OLS)
% Sunspots: u(t) = K + A sin(w t + phi) + e(t)
% Detrend:	y(t) = u(t) - E[u(t)]
% AR(2): 	y(t) = a y(t-1) + b y(t-2) + e(t)
Y = ssd(3:N);
Phi = [ssd(2:N-1), ssd(1:N-2)];         % Build the 'data matrix'
theta = pinv(Phi)*Y;                    % Estimate the parameters of the model
resid_var = var(Y - Phi*theta);         % Variance of the residuals
a = theta(1);
b = theta(2);

disp(sprintf('LS estimate of model parameters: a = %5.2f, b = %7.2f', a, b));
disp(sprintf('Variance of the residuals: %5.2f', resid_var));

%-------------------------------------------------------------------
% Estimate the power spectrum by PERIODOGRAM
dft = fft(ssd);                         % Discrete Fourier transform of the signal
dft = [ dft((N/2+1):N); dft(1:N/2); ];  % Interpret second chunk of the DFT as negative frequencies
dft_spectrum = abs(dft).^2./N;          % Periodogram

%-------------------------------------------------------------------
% Plot the spectrum estimates for COMPARISON
frequencies = linspace(-pi, pi, N);     % Frequency scale for plotting
ar_spectrum = resid_var./( 1 + a^2 + b^2 + 2*a*(b-1)*cos(frequencies) - 2*b*cos(2*frequencies) );	% Final step of a computation

figure('name', 'Spectrum Estimate'), title('Spectrum Estimate'), clf();
plot(frequencies, dft_spectrum, 'r'), hold on
plot(frequencies, ar_spectrum, 'k')
legend('Spectrum estimated with DFT', 'Spectrum estimated with AR fitting')
axis([-pi, pi]), xlabel('Frequency'), ylabel('Power spectral density'), hold off

%-------------------------------------------------------------------
% Compute the peak frequency and the phase of the (conjugate) poles of W(z),
% and the corresponding periods in years

freq_peak = acos(a*(b-1)/(4*b));
r = roots([1, -a, -b]);
freq_root = atan2( imag(r(1)), real(r(1)) );

printf('Peak Frequency:\t%5.2f\tCorresponding period (in years): %5.1f\n', freq_peak, 2*pi/freq_peak);
printf('Poles    Phase:\t%5.2f\tCorresponding period (in years): %5.1f\n', freq_root, 2*pi/freq_root);

Ysig = zeros(1:N, 1);
Ysig(1) = ss(1); Ysig(2) = ss(2); % suppose we got the first two datapoints.
for i=3:N
	Ysig(i) = a*Ysig(i-1)+b*Ysig(i-2)+ mean(ss)*rand(1);
end
figure(1), hold on, plot(yrs, Ysig, 'b')
legend('real data', 'AR fitting')

% order 2 is the minimum order such that the transfer function
% can have a pair of complex conjugate poles; a necessary condition for a resonance,
% that is a a “peak” in the power spectrum, and a more or less pronounced
% oscillatory behavior in realizations
