close all;
clear all;
rng(100);

dim = 7;
L = zeros(dim,1); U = zeros(dim,1);
dX = zeros(dim,1);

% Nominal values of the parameters
A = 7.049556277; B = 0.6022245584;
p = 4.0; q = 0.0; alpha = 1.80;
lambda = 21.0; gamma = 1.20;

N = [A;B;p;q;alpha;lambda;gamma];
L(:,1) = 0.9.*N(:); % lower-bound
U(:,1) = 1.1.*N(:); % upper-bound
ts = 30; % total number of samples
ns = 25; % subset of samples
U(4,1) = 0.1;

S = zeros(ts,dim);

for j = 1:dim
  S(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
  dX(j) = 1e-5.*(U(j,1) - L(j,1));
end

fid = fopen('params.txt','w');
fmt = '%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n';

for i = 1:ns
  fprintf(fid,fmt,S(i,1),S(i,2),S(i,3),S(i,4),S(i,5),S(i,6),S(i,7));
end

S_xpdx = zeros(ns,dim); 

for j = 1:dim
  S_xpdx = S(1:ns,:);
  S_xpdx(:,j) = S_xpdx(:,j) + dX(j);
  for i = 1:ns
    fprintf(fid,fmt,S_xpdx(i,1),S_xpdx(i,2),S_xpdx(i,3),S_xpdx(i,4),...
                    S_xpdx(i,5),S_xpdx(i,6),S_xpdx(i,7));
  end
end

fclose(fid);
