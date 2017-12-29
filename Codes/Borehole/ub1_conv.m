% This script computes and writes ubber bounds on total sobol indices for a
% range of sample size.

close all;
clear all;
rng(100);

fid = fopen('nsams_ub.txt','w');
fmt = '%d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n';
dim = 7;
dX = zeros(dim,1); pc = zeros(dim,1);
nu = zeros(dim,1); ub1 = zeros(dim,1);
L = [0.1;63070;990;63.1;700;1120;9855];
U = [0.2;115600;1110;116;820;1680;12045];

for i = 1:dim
    dX(i) = 1e-5.*(U(i)-L(i));
    pc(i) = ((U(i)-L(i)).^2)./(pi.^2);
end

pc(1) = 0.0161812.^2;

% Construct the set of samples based on the highest number i.e. 20
tsams = 500;
%S = zeros(tsams,dim);
% for j = 1:dim
%        S(:,j) = unifrnd(L(j,1),U(j,1),tsams,1);
% end

S = load('S.mat');

 % Evaluate model for the set S offline
G = uq_borehole(S.X);
 
for nsams = 5:5:tsams
    Gdx = zeros(nsams,dim);
    
    for i = 1:dim
        S_new = S.X(1:nsams,:);
        S_new(:,i) = S_new(:,i) + dX(i);
        Gdx(:,i) = uq_borehole(S_new);
        dGdx = (Gdx(:,i) - G(1:nsams,1))./dX(i);
        dGdx2 = dGdx.^2;
        T = sum(dGdx2);
        nu(i) = T./nsams;
        ub1(i) = pc(i).*nu(i);
    end
    ub1(1:dim) = ub1(1:dim)./sum(ub1);
    fprintf(fid,fmt,nsams,ub1(1),ub1(2),ub1(3),ub1(4),ub1(5),ub1(6),ub1(7));
    %fprintf(fid,'\n');
end

fclose(fid);
