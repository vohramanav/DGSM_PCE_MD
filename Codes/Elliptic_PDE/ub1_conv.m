% This script computes and writes ubber bounds on total sobol indices for a
% range of sample size.

close all;
clear all;
rng(100);

fid = fopen('nsams_ub.txt','w');
fmt = '%d %5.4f %5.4f %5.4f %5.4f\n';
dim = 4;
dX = zeros(dim,1); pc = zeros(dim,1);
nu = zeros(dim,1); ub1 = zeros(dim,1);
L = zeros(dim,1); U = zeros(dim,1);
L(1,1) = 0.09; U(1,1) = 0.11;
L(2:end,1) = 0.9; U(2:end,1) = 1.1;

for i = 1:dim
    dX(i) = 1e-5.*(U(i,1)-L(i,1));
    pc(i) = ((U(i,1)-L(i,1)).^2)./(pi.^2);
end

% Construct the set of samples based on the highest number i.e. 20
tsams = 20;
S = zeros(tsams,dim);
 for j = 1:dim
        S(:,j) = unifrnd(L(j,1),U(j,1),tsams,1);
 end

 % Evaluate model for the set S offline
 G = model_elliptic(S);
 
for nsams = 5:5:tsams
    Gdx = zeros(nsams,dim);
    
    for i = 1:dim
        S_new = S(1:nsams,:);
        S_new(:,i) = S_new(:,i) + dX(i);
        Gdx(:,i) = model_elliptic(S_new);
        dGdx = (Gdx(:,i) - G(1:nsams,1))./dX(i);
        dGdx2 = dGdx.^2;
        T = sum(dGdx2);
        nu(i) = T./nsams;
        ub1(i) = pc(i).*nu(i);
    end
    ub1(1:dim) = ub1(1:dim)./sum(ub1);
    fprintf(fid,fmt,nsams,ub1(1),ub1(2),ub1(3),ub1(4));
    %fprintf(fid,'\n');
end

fclose(fid);
