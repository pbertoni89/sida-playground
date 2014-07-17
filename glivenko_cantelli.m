%---------------------------------------------------------------------------------------
% Test of the Glivenko-Cantelli theorem
%---------------------------------------------------------------------------------------
% We run T times this experiment: draw a sequence of iid variables, plot the empirical distribution
% and the true distribution for comparison.
%
% (in abstract terms, we extract T "omega's", or "trajectories", and plot the empirical 
% distribution of their first N samples X1(omega) ... XN(omega) ).
%---------------------------------------------------------------------------------------

x = linspace(-5, 5, 2000);        % Abscissas for plotting
trueDist = normcdf(x);            % Cumulative distribution of N(0,1)

T = 20;                            % Number of trajectories
MAX_N = 10000;                    % Length of the trajectory (ideally infinite)
trajectories = randn(T, MAX_N);   % (Ideally) iid variables distributed as N(0,1)
i = 1;

% Show in subsequent plots what happens as N increases, to appreciate uniform convergence
for N=[5, 10, 30, 100, 300, 1000, 3000, 10000],

   figure(i);
   clf;
   hold on;
   axis([-5 5 0 1]);
   title(sprintf('N = %d', N));

   % For each trajectory, build and plot the graph of the empirical distribution up to time N
   for k = 1:T,
      empiricalDist = zeros(1, length(x));
      for i = 1:N,
         stepFunction = trajectories(k, i) <= x;          % The indicator function 1_{Xi <= x}
         empiricalDist = empiricalDist + stepFunction/N;  % Accumulates weighted sum
      end
      plot(x, empiricalDist, 'r');
   end

   % Then plot the true one
   plot(x, trueDist, 'k');
end

%---------------------------------------------------------------------------------------
% Now try with T = 2, 10, 100 trajectories.
% The empirical distribution converges UNIFORMLY to the true one
% for "almost all" trajectories (i.e. ALMOST SURELY, in the abstract setting,
% but in practice for all of them!).
%---------------------------------------------------------------------------------------
