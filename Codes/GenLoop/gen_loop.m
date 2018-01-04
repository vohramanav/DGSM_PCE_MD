function [t y] = gen_loop(Tf)
% Solves the genetic positive feed-back loop system

% reaction rates
c = [50 1000 50 1000 1 10 3 1 6];

% solve system
[t y] = gen_loop_solver(c, Tf);

% plotting
figure(1)
plot(t, y(:,2), 'linewidth', 2);




% ------------------------------------------
% Sub-function for solving the RRE system
function [t y] = gen_loop_solver(c, Tf)

% give the reaction rates their physical names
kappa_plus  = c(1)/2;
kappa_minus = c(2); 
k_plus      = c(3);
k_minus     = c(4);
alpha       = c(5); 
beta        = c(6);
sigma       = c(7); 
gamma_p     = c(8); 
gamma_m     = c(9);

% Note: Y(1) = x, 
%       Y(2) = y, 
%       Y(3) = d0, 
%       Y(4) = dr, 
%       Y(5) = m

% ------------------------------------
% define the right hand side of system
% ------------------------------------
xdot  = @(Y)(2*kappa_minus*Y(2) - 2*kappa_plus*Y(1)^2 + sigma*Y(5) - gamma_p*Y(1));
ydot  = @(Y)(kappa_plus * Y(1)^2 - kappa_minus*Y(2) + k_minus*Y(4) - k_plus*Y(3)*Y(2));
d0dot = @(Y)(k_minus*Y(4) - k_plus*Y(2)*Y(3));
drdot = @(Y)(k_plus*Y(2)*Y(3) - k_minus*Y(4));
mdot  = @(Y)(alpha*Y(3) + beta*Y(4) - gamma_m*Y(5));

rhs = @(t, Y)([xdot(Y);
               ydot(Y);
               d0dot(Y); 
               drdot(Y);
               mdot(Y)]);

% solve the ODE:
Y0 = [10 0 20 0 0];
%tmesh = [0 : 1e-6 : 0.2, 0.3 : 0.01 : Tf];
tmesh = [0 Tf];
[t, y] = ode15s(rhs, tmesh, Y0);


