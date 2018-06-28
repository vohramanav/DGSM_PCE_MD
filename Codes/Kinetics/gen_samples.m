close all;
clear all;
rng(40);

dim = 19;
L = zeros(1,dim); U = zeros(1,dim); dX = zeros(1,dim);

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);
nsams = 10; pts = zeros(nsams*(dim+1),dim);

for j = 1:dim
  pts(1:nsams,j) = unifrnd(L(1,j),U(1,j),nsams,1);
  dX(1,j) = 1e-3.*(U(1,j) - L(1,j));
end

for i = 1:dim
  pts(i*nsams+1:(i+1)*nsams,:) = pts(1:nsams,:);
  pts(i*nsams+1:(i+1)*nsams,i) = pts(i*nsams+1:(i+1)*nsams,i) + dX(1,i);
end

f1 = fopen('pts_test.txt','w');
%fmt = '%15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f\n';
fmt = '%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n';

for i = 1:nsams*(dim+1)
  fprintf(f1,fmt,pts(i,1),pts(i,2),pts(i,3),pts(i,4),pts(i,5),pts(i,6),pts(i,7),pts(i,8),pts(i,9),pts(i,10),pts(i,11),pts(i,12),pts(i,13),pts(i,14),pts(i,15),pts(i,16),pts(i,17),pts(i,18),pts(i,19));
end

fclose(f1);
