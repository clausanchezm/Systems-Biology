function dydt = lotka_volterra(t,y,alpha,beta)

% Lotka-Volterra predator-prey model: a system of coupled differential
% equations. 
% This is a function that returns a column vector of state derivatives,
% given state, time, and parameter values.

dydt = zeros(2,1);

dydt(1) = (1 - alpha*y(2))*y(1);
dydt(2) = (-1 + beta*y(1))*y(2);

end