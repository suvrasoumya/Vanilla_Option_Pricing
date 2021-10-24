%%
clc
clear all
close all
%%
call = CrankNicGauss(100,100,0.03, 0.3, 15, 0.01, 1)
% CrankNicGauss(Strike, Stock price, interest_rate, volatility, delta_S, delta_t, T_maturity)
[Call_bsPrice, not_used] = blsprice(100, 100, 0.03, 1, 0.3, 0);
%%
function oPrice = CrankNicGauss(X,S0,r,sig, del_S,del_t, time)
% Function to calculate the price of a vanilla European

% Create mesh 

S_vector = 0:del_S:3*S0;
t_vector = 0:del_t:time;

% Specify the boundary conditions for call

% Get the number of grid points
M = length(S_vector)-1;
N = length(t_vector)-1;

% Pre-allocate the output
price_mesh(1:M+1,1:N+1) = nan;

price_mesh(:,end) = max(S_vector-X,0);
price_mesh(1,:) = 0;
price_mesh(end,:) = (S_vector(end)-X)*exp(-r*t_vector(end:-1:1));

% Calculate the coefficients of tridiagonal matrix
j = 0:M;
sig2 = sig*sig;
pj = (del_t/4)*(sig2*(j.^2) - r*j);
qj = -(del_t/2)*(sig2*(j.^2) + r);
rj = (del_t/4)*(sig2*(j.^2) + r*j);

% The tridiagonal matrix for LHS and RHS are C and D respectively
C = -diag(pj(3:M),-1) + diag(1-qj(2:M)) - diag(rj(2:M-1),1);
D = diag(pj(3:M),-1) + diag(1+qj(2:M)) + diag(rj(2:M-1),1);

% Solve at each node
% bound_S function is used to fetch the boundary condition for the stock
% at +/- infinity for option price

bound_S = zeros(size(D,2),1);
dx = zeros(size(D,2),1);
for idx = N:-1:1
    if length(bound_S)==1
        bound_S = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1)) + ...
            rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    else
        bound_S(1) = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1));
        bound_S(end) = rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    end
    price_mesh(2:M,idx) = gauss_seidel(C,(D*price_mesh(2:M,idx+1) + bound_S), 200, 1.e-6);
end

% Calculate the option price
oPrice = interp1(S_vector,price_mesh(:,1),S0);
end

% Reference Taking Gauss-Siedel pseudo-code from University of Waterloo
% Publication.
function x = gauss_seidel( M, b, N, e )
	% Solve Mx = b
	% The diagonal entries of M and their inverses
	n = length( b );
	d = diag( M );

	if ~all( d ) 
		error 'at least one diagonal entry is zero';
	end

	invd = d.^-1;
	% Matrix of off-diagonal entires of N
	Moff = M - diag( d );

	% Use d.^-1*b as the first approximation to x
	invdb = invd.*b;
	x = d.*b;

	%              -1
	% Iterate x = D  (b - M   *x)
	%                      off
	for k = 1:N
		xprev = x;

		for i = 1:n
			x(i) = invdb(i) - invd(i).*(Moff(i,:)*x);
		end

		if norm( x - xprev, inf ) < e
			return;
		end
	end

	error 'the method did not converge';
end