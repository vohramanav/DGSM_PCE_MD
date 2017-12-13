function Umax = semilinear_elliptic_pde(k,c,al,be)
% Solve the problem
%
%     -kappa Laplacian u + c u^3 = q
%     + zero Dirichlet BCs
%       
%     on [0,1] x [0,1].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using Newton's method.
%
%     we use kappa = .1, and q(x) defined by
%
%
%     q(x,y) = (-2pi^2)(cos(2 pi x)sin(pi y)^2 + sin(pi x)^2 * cos(2 pi y))
%
% the mesh resolution level (m) and coefficient in front of the nonlinear term
% (c) are input.


a = 0; 
b = 1; 
m = 1e3;
h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(a,b,m+2);   % grid points y including boundaries


[X,Y] = ndgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:m+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);


frhs = @(x,y)((-2*pi^2)*(al.*cos(2*pi*x) .* sin(pi*y).^2 + be.*sin(pi*x).^2 .* cos(2*pi*y)));


rhs = frhs(Xint,Yint);        % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.

% utrue = sin(pi*x).^2 .* sin(pi*y).^2;

usoln = zeros(size(X));              % use true solution for this test problem
                            % This sets full array, but only boundary values
                            % are used below.  For a problem where utrue
                            % is not known, would have to set each edge of
                            % usoln to the desired Dirichlet boundary values.


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
rhs(:,m) = rhs(:,m) - usoln(Iint,m+2)/h^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;

% convert the 2d grid function rhs into a column vector for rhs of system:
q = reshape(rhs,m*m,1);

% form matrix A:
%kappa = .1;
kappa = k;
I = speye(m);
e = ones(m,1);
T = spdiags([e -4*e e],[-1 0 1],m,m);
S = spdiags([e e],[-1 1],m,m);
A = -kappa*(kron(I,T) + kron(S,I)) / h^2;

%
% setup and solve the nonlinear system
%
F = @(u)(A * u + c * u.^3 - q);
DF = @(u)(A + spdiags(3*c*u.^2, 0, m^2, m^2)); 
%DF = @(u)(A + diag(3*c*u.^2));


% tolerances
atol =1e-10; rtol=1e-10; 

% iteration data
maxit = 20;
it = 0; 

% initial guess 
u0 = zeros(m^2,1);
res_norm0 = norm(F(u0));
res_norm = res_norm0;
resvec = res_norm0;

u = u0;
while (res_norm > (atol + rtol*res_norm0)) && (it < maxit)

   Fxc = F(u);
   DFxc = DF(u);
   step = -DFxc \ Fxc; 
   u = u + step; 

   % record iteration history 
   Fxp = F(u);
   res_norm = norm(Fxp); 
   resvec = [resvec,res_norm]; 
  
   it = it + 1;
end

Umax = max(u);
save('umax.mat','Umax');
usoln(Iint,Jint) = reshape(u,m,m);

% plot results:
close all;

% figure
% semilogy(resvec,'-o', 'linewidth', 2);
% set(gca, 'fontsize', 20);
% xlabel('k');
% ylabel('relative residual');
 
 % plot solution:
%  figure
%  surf(X,Y,usoln);
%  xlabel('$$\mathrm{x}$$','interpreter','latex');
%  ylabel('$$\mathrm{y}$$','interpreter','latex');
%  set(gca,'fontsize',14);
%  set(gca,'TickLabelInterpreter','latex');
%  set(gcf,'color',[1,1,1]);
%  shading interp;
%  axis square;
%  box on;grid on;
%  c = colorbar();
%  set(c,'TickLabelInterpreter','latex');
%  title('$$\mathrm{\kappa = 0.1~~c = 1.0~~a = 1.0~~b = 1.0}$$','interpreter','latex')
%  print -dpng u_soln.png 
end
