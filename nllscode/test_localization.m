% Problem: localization of a transmitter with nonlinear least-squares

function test_localization();

	global TX_POWER;
	TX_POWER = 30;         % Transmission power
	ERROR_VARIANCE = 0.01; % Variance of power measurement errors at each receiver

	% Positions of the receiving stations on the XY plane
	stations = [
	              0,  0;	%
	              0, 10;	%
	             10,  0;	%
	             10, 10; 	%
	              5,  5;	%
	             15, 15;
	           ];
	           
	% True position of the transmitter on the XY plane (= theta_zero)
	transmitter_unknown = [6; 6];  
	
	N = size(stations, 1);
	Y = zeros(N,1);
	for i=1:N,
		% Measures corrupted by noise
		Y(i) = received_power(stations(i,:)', transmitter_unknown) + sqrt(ERROR_VARIANCE)*randn(1,1);
	end

	% Estimate the position of the transmitter with NLLS
	initial_guess = [1;1];
	transmitter_est = levenberg_marquardt(Y, stations, @received_power, @received_power_grad, initial_guess)

	figure(1), clf(), scatter(stations(:,1), stations(:,2), 40, 'r', "filled")
	hold on
	scatter(transmitter_unknown(1), transmitter_unknown(2), 40, 'c', "filled")
	scatter(transmitter_est(1), transmitter_est(2), 30, 'g', "filled")

end;

% (NON LINEAR) Function: power received from a transmitter (position = theta = (x,y)) 
% by a station (position = (x_i,y_i))
% f(x,y) = TX_POWER / ( (x_i - x)^2 + (y_i - y)^2 )
function pow = received_power(station_pos, transmitter_pos);
	global TX_POWER;
	distance = norm(station_pos - transmitter_pos);
	pow = TX_POWER/distance^2;
end;

% Compute the derivative with respect to theta = (x,y) of the function
% f(x,y) = TX_POWER / ( (x_i - x)^2 + (y_i - y)^2 )
function grad = received_power_grad(station_pos, transmitter_pos);
	global TX_POWER;
	distance = norm(station_pos - transmitter_pos);

	df_dx = 2*TX_POWER*(station_pos(1) - transmitter_pos(1)) / distance^4;
	df_dy = 2*TX_POWER*(station_pos(2) - transmitter_pos(2)) / distance^4;
	grad = [df_dx, df_dy];
end;
