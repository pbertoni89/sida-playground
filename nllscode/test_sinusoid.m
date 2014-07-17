% Problem: find amplitude, frequency and phase of a sinusoid
% with nonlinear least squares
% Now try with different initial guesses.
% Note how the estimated model is particularly influenced
% by the guess of amplitude and frequency

function test_sinusoid();

	N = 100; % Number of measures
	ERROR_VARIANCE = 0.01; % Variance of measurement errors
	close all

	%---------------------------------------------
	% Generate random measures	
	 T = linspace(0, 1, N)';
	%T = rand(N, 1);                % Random times
	theta_true = [1; 3.5; pi/4];   % True parameter, theta_zero = (A, F, phi)

	Y = zeros(N,1);
	for i=1:N,
		% Measures corrupted by noise
		Y(i) = sinusoid(T(i), theta_true) + sqrt(ERROR_VARIANCE)*rand(1,1);
	end

	%---------------------------------------------
	% Estimate amplitude, frequency and phase
	GUESS = 5;
	%initial_guess(GUESS, 3);
	for j=1:GUESS
		initial_guess = rand(3,1)*4; %[1.3; 4; 0];
		theta_est = levenberg_marquardt(Y, T, @sinusoid, @sinusoid_grad, initial_guess);
		% Note: the phase is meant to be correct only modulo 2*pi

		Nplot = 200;
		t = linspace(0, 1, Nplot);
		yguess = zeros(1, Nplot);
		for i=1:Nplot,
			yguess(i) = sinusoid(t(i), initial_guess);
		end

		ypred = zeros(1, Nplot);
		for i=1:Nplot,
			ypred(i) = sinusoid(t(i), theta_est);
		end
		
		figure(j), clf, hold on
		plot(T, Y', 'b*'), plot(t, yguess, 'r'), plot(t, ypred, 'k')
		axis([0 1 -4 4]), xlabel('Time'), ylabel('Value');
		legend('Measures', 'Initial guess', 'Estimated model')
		title('Estimation of a sinusoidal model')
		
		printf("guessing <%f,%f,%f> from initial <%f,%f,%f> produces <%f,%f,%f>\n", 
			theta_true(1), theta_true(2), theta_true(3), 
			initial_guess(1), initial_guess(2), initial_guess(3), 
			theta_est(1), theta_est(2), theta_est(3))
	end
end;

% Function: y = A*sin(2*pi*F*t + phi);
function y = sinusoid(t, params);
	A 	= params(1);
	F 	= params(2);
	phi = params(3);
	y 	= A*sin(2*pi*F*t + phi);
end;

% Compute the derivative with respect to theta = (A, F, phi) of the function
function grad = sinusoid_grad(t, params);
	A 	= params(1);
	F 	= params(2);
	phi = params(3);

	df_dA 	= 	sin(2*pi*F*t + phi);
	df_dF 	= A*cos(2*pi*F*t + phi)*2*pi*t;
	df_dphi = A*cos(2*pi*F*t + phi);

	grad = [df_dA, df_dF, df_dphi];
end;
