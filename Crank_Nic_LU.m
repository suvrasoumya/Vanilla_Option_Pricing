%%
clc
clear all
close all
%%
call = CrankNicLU(60, 100, 0.03, 0.3, 2, 0.01, 1)
% CrankNicGauss(Strike, Stock price, volatility, delta_S, delta_t, T_maturity)
[Call_bsPrice, not_used] = blsprice(100,60, 0.05, 1, 0.2, 0);
%%
function call_price = CrankNicLU(X, S0, r, sig, del_S, del_t, time)

% Create mesh 

S_vector = 0:del_S:2*S0;
t_vector = 0:del_t:time;

% Get the number of grid points
M = length(S_vector)-1;
N = length(t_vector)-1;

% Pre-allocate the output
price_mesh(1:M+1,1:N+1) = nan;

% Specify the boundary conditions
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
[L,U] = lu(C);
D = diag(pj(3:M),-1) + diag(1+qj(2:M)) + diag(rj(2:M-1),1);

% Solve at each node
% bound_S function is used to fetch the boundary condition for the stock
% at +/- infinity for option price

bound_S = zeros(size(D,2),1);
for idx = N:-1:1
    if length(bound_S)==1
        bound_S = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1)) + ...
            rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    else
        bound_S(1) = pj(2)*(price_mesh(1,idx)+price_mesh(1,idx+1));
        bound_S(end) = rj(end)*(price_mesh(end,idx)+price_mesh(end,idx+1));
    end
    price_mesh(2:M,idx) = U\(L\(D*price_mesh(2:M,idx+1) + bound_S));
end

% Calculate the option price
call_price = interp1(S_vector,price_mesh(:,1),S0);
end