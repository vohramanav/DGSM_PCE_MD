close all;
clear all;
rng(40);

dim = 5;
L = zeros(dim,1); U = zeros(dim,1);

% Nominal values of the parameters
A = 7.049556277; B = 0.6022245584;
p = 4.0; q = 0.0; alpha = 1.80;
lambda = 21.0; gamma = 1.20;

N = [A;q;alpha;lambda;gamma];
L(:,1) = 0.9.*N(:); % lower-bound
U(:,1) = 1.1.*N(:); % upper-bound
ts = 1e6; % total number of samples
U(2,1) = 0.1;

S = zeros(ts,dim);

for j = 1:dim
  S(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
end

fid = fopen('params.txt','w');
fmt = '%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n';

for i = 1:ts
  fprintf(fid,fmt,S(i,1),B,p,S(i,2),S(i,3),S(i,4),S(i,5));
end

fclose(fid);
