% Model: y(t) = a y(t-1) + b u(t-1) + e(t)
% e(t) is a white noise; u(t) and y(t) are measured without errors

nruns = 100;
a_estimates = zeros(1, nruns);
b_estimates = zeros(1, nruns);

for run=1:nruns,
	% Construction of the process
	N = 1000;                   % Time horizon
	a = 0.8;                    % Parameter
	b = 0.2;                    % Parameter
	u = randn(N,1);             % Some input signal with persistent excitation
	e = 0.3*randn(N+1,1);       % Process noise
	y = zeros(N+1, 1);
	y(1) = 4;                   % Initial condition
	for t=2:N+1,
		y(t) = a*y(t-1) + b*u(t-1) + e(t);
	end

	% Data for least squares
	ypast    = y(1:N);			% 1234567....1000
	ypresent = y(2:N+1);		%  234567........1001
  	Phi = [ypast, u];           % Regressors. size(Phi)= {N, 2}

	% Least squares estimation
	theta_LS = pinv(Phi)*ypresent;    % Pseudoinverse. Same as: inv(Phi'*Phi)*Phi'*ypresent
	a_estimates(run) = theta_LS(1);
	b_estimates(run) = theta_LS(2);
	
	%Phi3 = [ypast.^2, ypast, u];
	%theta_LS3 = pinv(Phi3)*ypresent;
	%a_estimates3(run) = theta_LS3(1);
	%b_estimates3(run) = theta_LS3(2);
	%c_estimates3(run) = theta_LS3(3);
end

disp(sprintf('LS estimate of a over %d runs: average %7.5f, variance %7.5f', 
             nruns, mean(a_estimates), var(a_estimates)));
disp(sprintf('LS estimate of b over %d runs: average %7.5f, variance %7.5f\n',
             nruns, mean(b_estimates), var(b_estimates)));

%disp(sprintf('LS3 estimate of a over %d runs: average %7.5f, variance %7.5f', 
%             nruns, mean(a_estimates3), var(a_estimates3)));
%disp(sprintf('LS3 estimate of b over %d runs: average %7.5f, variance %7.5f',
%             nruns, mean(b_estimates3), var(b_estimates3)));
%disp(sprintf('LS3 estimate of c over %d runs: average %7.5f, variance %7.5f', 
%             nruns, mean(c_estimates3), var(c_estimates3)));

% l'idea è, per ogni istante del PS, i.e. per ogni VC, effettuare un'identificazione OLS.
% successivamente stampare la media e la varianza temporale per ogni parametro stimato.
