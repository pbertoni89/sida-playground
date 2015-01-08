%---------------------------------------------------------------------------------------
% Test of the Glivenko-Cantelli theorem
%---------------------------------------------------------------------------------------
% We run T times this experiment: draw a sequence of iid variables, plot the empirical distribution
% and the true distribution for comparison.
%
% (in abstract terms, we extract T "omega's", or "trajectories", and plot the empirical 
% distribution of their first N samples X1(omega) ... XN(omega) ).
%---------------------------------------------------------------------------------------
clc, close all

X_LEN = 2000;
X_EDGE = 5;
x = linspace(-X_EDGE, X_EDGE, X_LEN);        % Abscissas for plotting
trueDist = normcdf(x);            % Cumulative distribution of N(0,1) in space x

T = 10;                            % Number of trajectories
MAX_N = 10000;                    % Length of the trajectory (ideally infinite)
trajectories = randn(T, MAX_N);   % (Ideally) iid variables distributed as N(0,1)

% Show in subsequent plots what happens as N increases, to appreciate uniform convergence
for N=[5, 10, 30, 100, 300, 1000, 3000, 10000],

   figure(), clf, hold on, axis([-X_EDGE X_EDGE 0 1]), title(sprintf('N = %d', N)), xlabel('x'), ylabel('F(x)')

   % For each trajectory, build and plot the graph of the empirical distribution up to time N
   for k = 1:T,
      empiricalDist = zeros(1, X_LEN);
      for i = 1:N,
         stepFunction = trajectories(k, i) <= x;         % The indicator function 1_{Xi <= x}
         empiricalDist = empiricalDist + stepFunction/N;  % Accumulates weighted sum
      end
      plot(x, empiricalDist, 'r')
   end
   plot(x, trueDist, 'k')	% Then plot the true one
end

%---------------------------------------------------------------------------------------
% Now try with T = 2, 10, 100 trajectories.
% The empirical distribution converges UNIFORMLY to the true one
% for "almost all" trajectories (i.e. ALMOST SURELY, in the abstract setting,
% but in practice for all of them!).
%---------------------------------------------------------------------------------------
