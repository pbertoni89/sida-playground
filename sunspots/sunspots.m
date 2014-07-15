% Octave/Matlab code:
% Autoregressive fit of Wolf's sunspot numbers
% Comparison of the power spectral densities estimates obtained from DFT and AR fit 

%-------------------------------------------------------------------
% Load data
sunspot_data = load('sunspots.txt');
N = size(sunspot_data,1);
yrs = sunspot_data(:, 1);     % Years
ss  = sunspot_data(:, 2);     % Wolf's sunspot numbers
ssd = ss - mean(ss);          % Detrend the time series

%-------------------------------------------------------------------
% Plot sunspot numbers
figure(1);
clf();
plot(yrs, ss, 'r');
axis([1749, 1927, 0, 160]);
xlabel('year');
ylabel('Wolf''s number');
%print -color -depslatexstandalone wolf.eps

%-------------------------------------------------------------------
% Estimate the power spectrum by PEM (ordinary least squares)
% Model: y(t) = a y(t-1) + b y(t-2) + e(t)
Y = ssd(3:N);
Phi = [ssd(2:N-1), ssd(1:N-2)];         % Build the 'data matrix'
theta = pinv(Phi)*Y;                    % Estimate the parameters of the model
resid_var = var(Y - Phi*theta);         % Variance of the residuals
a = theta(1);
b = theta(2);

disp(sprintf('LS estimate of model parameters: a = %5.2f, b = %7.2f', a, b));
disp(sprintf('Variance of the residuals: %5.2f', resid_var));

%-------------------------------------------------------------------
% Estimate the power spectrum with the periodogram
dft = fft(ssd);                         % Discrete Fourier transform of the signal
dft = [ dft((N/2+1):N); dft(1:N/2); ];  % Interpret second chunk of the DFT as negative frequencies
dft_spectrum = abs(dft).^2./N;          % Periodogram

%-------------------------------------------------------------------
% Plot the spectrum estimates for comparison
frequencies = linspace(-pi, pi, N);     % Frequency scale for plotting
ar_spectrum = resid_var./( 1 + a^2 + b^2 + 2*a*(b-1)*cos(frequencies) - 2*b*cos(2*frequencies) );

figure(2);
clf();
plot(frequencies, dft_spectrum, 'r');
hold on;
plot(frequencies, ar_spectrum, 'k');
legend('Spectrum estimated with DFT', 'Spectrum estimated with AR fitting');
axis([-pi, pi]);
xlabel('Frequency');
ylabel('Power spectral density');
hold off;
%print -color -depslatexstandalone wolfspec.eps

%-------------------------------------------------------------------
% Compute the peak frequency and the phase of the (conjugate) poles of W(z),
% and the corresponding periods in years

freq_peak = acos(a*(b-1)/(4*b));
r = roots([1, -a, -b]);
freq_root = atan2( imag(r(1)), real(r(1)) );

disp(sprintf('Peak frequency: %5.2f', freq_peak));
disp(sprintf('Corresponding period (in years): %5.1f', 2*pi/freq_peak));
disp(sprintf('Phase of the poles: %5.2f', freq_root));
disp(sprintf('Corresponding period (in years): %5.1f', 2*pi/freq_root));



