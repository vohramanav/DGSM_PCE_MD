close all;
clear all;
rng(100);

dim = 7;
dX = zeros(dim,1); pc = zeros(dim,1);
mu = zeros(dim,1); ub = zeros(dim,1);
L = zeros(dim,1); U = zeros(dim,1);

% Nominal values of the parameters
A = 7.049556277; B = 0.6022245584;
p = 4.0; q = 0.0; alpha = 1.80;
lambda = 21.0; gamma = 1.20;

N = [A;B;p;q;alpha;lambda;gamma];
L(:,1) = 0.9.*N(:); % lower-bound
U(:,1) = 1.1.*N(:); % upper-bound

for j = 1:dim
  dX(i) = 1e-5.*(U(i,1)-L(i,1));
  pc(j) = ((U(i,1)-L(i,1)).^2)./(pi.^2);
end

nsams = 5;
G = zeros(nsams,1); Gdx = zeros(nsams,dim);
G(:,1) = k_data(1:nsams,2);

for j = 1:dim
  Gdx(:,j) = k_data(j*nsams+1:(j+1)*nsams,2);
  dGdx = (Gdx(:,j) - G(:,1))./dX(j);
  dGdx2 = dGdx.^2;
  T = sum(dGdx2);
  nu(j) = T./nsams;
  ub(j) = pc(j).*nu(j);
end

ub(1:dim) = ub(1:dim)./sum(ub);

