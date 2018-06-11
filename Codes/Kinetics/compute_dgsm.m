close all;
clear all;

dim = 19;
L = zeros(1,dim); U = zeros(1,dim); dX = zeros(1,dim);
pc = zeros(1,dim); mu = zeros(1,dim); ub = zeros(1,dim);
id_data_rng10 = load('record_id_rng10.txt');
id_data_rng20 = load('record_id_rng20.txt');
id_data_rng30 = load('record_id_rng30.txt');
id_data_rng40 = load('record_id_rng40.txt');

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);

for j = 1:dim
  dX(1,j) = 1e-3.*(U(1,j) - L(1,j));
  pc(1,j) = ((U(1,j)-L(1,j)).^2)./(pi.^2);
end

nsams = 20;
G = zeros(nsams,1); Gdx = zeros(nsams,dim);
G(1:5,1) = id_data_rng10(1:5,1);
G(6:10,1) = id_data_rng20(1:5,1);
G(11:15,1) = id_data_rng30(1:5,1);
G(16:20,1) = id_data_rng40(1:5,1);
%save('id_7D_20pts.dat','G','-ASCII')

for j = 1:dim
  Gdx(1:5,j) = id_data_rng10(j*5+1:(j+1)*5,1);
  Gdx(6:10,j) = id_data_rng20(j*5+1:(j+1)*5,1);
  Gdx(11:15,j) = id_data_rng30(j*5+1:(j+1)*5,1);
  Gdx(16:20,j) = id_data_rng40(j*5+1:(j+1)*5,1);
  dGdx = (Gdx(:,j) - G(:,1))./dX(1,j);
  dGdx2 = dGdx.^2;
  T = sum(dGdx2);
  mu(1,j) = T./nsams;
  ub(1,j) = pc(1,j).*mu(1,j);
end

ub(1,1:dim) = ub(1,1:dim)./sum(ub);
save('ub20.mat','ub');
save('mu20.mat','mu');
