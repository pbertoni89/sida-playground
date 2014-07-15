% Model: y(t) = a y(t-1) + b u(t-1)
% u(t) is measured without errors
% y(t) is measured with an error, ym(t) = y(t) + e(t)
 
nruns = 100;   % Number of identification experiments
a_estimates = zeros(1, nruns);
b_estimates = zeros(1, nruns);

for run=1:nruns,
   % Construction of the process
   N = 1000;                   % Time horizon
   a = 0.8;                    % Parameter
   b = 0.2;                    % Parameter
   u = randn(N,1);             % Some input signal with persistent excitation
   e = 0.3*randn(N+1,1);       % Some random noise
   y = zeros(N+1, 1);
   y(1) = 4;                   % Initial condition
   for t=2:N+1,
      y(t) = a*y(t-1) + b*u(t-1);
   end
   ym = y + e;                 % Measured output

   % Data for least squares
   ypast    = ym(1:N);
   ypresent = ym(2:N+1);
   Phi = [ypast, u];           % Regressors

   % Least squares estimation
   theta_LS = pinv(Phi)*ypresent;
	a_estimates(run) = theta_LS(1);
	b_estimates(run) = theta_LS(2);
end

disp(sprintf('LS estimate of a over %d runs: average %7.5f, variance %7.5f',
             nruns, mean(a_estimates), var(a_estimates)));
disp(sprintf('LS estimate of b over %d runs: average %7.5f, variance %7.5f',
             nruns, mean(b_estimates), var(b_estimates)));

