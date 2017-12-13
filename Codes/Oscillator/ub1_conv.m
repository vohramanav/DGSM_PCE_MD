% This script computes and writes ubber bounds on total sobol indices for a
% range of sample size.

close all;
clear all;
rng(100);

mu = [1.0,1.0,0.1,0.5,1.0,1.0];
sigma = [0.05.^2,0,0,0,0,0;...
    0,0.1.^2,0,0,0,0;...
    0,0,0.01.^2,0,0,0;...
    0,0,0,0.05.^2,0,0;...
    0,0,0,0,0.2.^2,0;...
    0,0,0,0,0,0.2.^2];

fid = fopen('nsams_ub.txt','w');
fmt = '%d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n';
dim = 6;
dX = zeros(dim,1); pc = zeros(dim,1);
nu = zeros(dim,1); ub1 = zeros(dim,1);

for i = 1:dim
    dX(i) = 1e-5.*mu(i);
    pc(i) = sigma(i,i);
end

for nsams = 5:5:500
    X = lhsnorm(mu,sigma,nsams);
    G = model_oscillator(X);
    Gdx = zeros(nsams,dim);
    
    for i = 1:dim
        X_new = X;
        X_new(:,i) = X_new(:,i) + dX(i);
        Gdx(:,i) = model_oscillator(X_new);
        dGdx = (Gdx(:,i) - G)./dX(i);
        dGdx2 = dGdx.^2;
        S = sum(dGdx2);
        nu(i) = S./nsams;
        ub1(i) = pc(i).*nu(i);
    end
    ub1(1:dim) = ub1(1:dim)./sum(ub1);
    fprintf(fid,fmt,nsams,ub1(1),ub1(2),ub1(3),ub1(4),ub1(5),ub1(6));
    %fprintf(fid,'\n');
end

fclose(fid);
