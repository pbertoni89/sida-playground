%--------------------------------------------------------------
% Levenberg-Marquardt method
% measurement model: y = f(x, theta) + eps
%--------------------------------------------------------------
% Input:
% Y: measures from 1 to N, stacked in a COLUMN
% X: explanatory data from 1 to N, stacked vertically (one x_i per each row)
%
% F_function: the function f(x, theta).
% Must accept both x and theta as column vectors and return a scalar
%
% F_gradient: the gradient of f(x, theta).
% Must accept both x and theta as column vectors and return the gradient w.r.t. theta as a ROW
%--------------------------------------------------------------

function thetaopt = levenberg_marquardt(Y, X, F_function, F_gradient, theta0, TOL = 1e-10);
	TAU = 1e-3;
	MAX_ITERATIONS = 200;
	iterations = 0;

	p = length(theta0);
	m = size(X, 2);
	N = size(X, 1);
	assert(N == length(Y));
	theta = theta0;

	% Initialize lambda
	nu = 2;
	lambda = TAU*max(diag(J'*J));

	loop = 1;
	while loop,

		% Build current F(theta) and J(theta)
		F = build_F_vector(F_function, X, theta);
		J = build_F_jacobian(F_gradient, X, theta);

		currentQ = norm(Y - F)^2;
		JtY_F = J'*(Y - F);

		reject_displacement = 1;
		while reject_displacement,

			% Compute Levenberg-Marquadt displacement
			dtheta = pinv(J'*J + lambda*eye(p))*JtY_F;

			% Compute F(theta + dtheta)
			newF = build_F_vector(F_function, X, theta + dtheta);
			newQ = norm(Y - newF)^2;

			predicted_decrease = dtheta'*(JtY_F + lambda*dtheta);
			if predicted_decrease > 0,
				% Gain ratio
				rho = (currentQ - newQ) / predicted_decrease;

				if rho > 0,
					% accept the step
					theta = theta + dtheta;
					reject_displacement = 0;
					lambda = lambda*max(1/3, 1 - (2*rho - 1)^3);
					nu = 2;
				else
					lambda = nu*lambda;
					nu = 2*nu;
				end
			else
				reject_displacement = 0;
			end
		end

		% Check: current solution is acceptable or we have entered an
		% endless loop (pathological condition)
		iterations = iterations + 1;
		if norm(dtheta) <= TOL || iterations >= MAX_ITERATIONS,
			loop = 0;
		end
	end

	thetaopt = theta;
end


function F = build_F_vector(F_function, X, theta);
	N = size(X, 1);
	F = zeros(N, 1);
	for i = 1:N,
		F(i) = F_function(X(i,:)', theta);
	end
end

function J = build_F_jacobian(F_gradient, X, theta);
	N = size(X, 1);
	p = length(theta);
	J = zeros(N, p);
	for i = 1:N,
		J(i,:) = F_gradient(X(i,:)', theta);
	end
end


