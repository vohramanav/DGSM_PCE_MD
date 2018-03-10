close all;
clear all;

k1 = load('ref_k/k_1_80.txt');
k2 = load('ref_k/k_81_160.txt');
k3 = load('ref_k/k_161_240.txt');
p1 = load('ref_params/ref_p_1_80.txt');
p2 = load('ref_params/ref_p_81_160.txt');
p3 = load('ref_params/ref_p_161_240.txt');
sk1 = size(k1,1);
sk2 = size(k2,1);
sk3 = size(k3,1);
sp1 = size(p1,1);
sp2 = size(p2,1);
sp3 = size(p3,1);
kc = zeros(sk1+sk2+sk3,1);
pc = zeros(sp1+sp2+sp3,7);
kc(1:sk1,1) = k1(:,2);
kc(sk1+1:sk1+sk2,1) = k2(:,2);
kc(sk1+sk2+1:sk1+sk2+sk3,1) = k3(:,2);
pc(1:sp1,:) = p1(:,:);
pc(sp1+1:sp1+sp2,:) = p2(:,:);
pc(sp1+sp2+1:sp1+sp2+sp3,:) = p3(:,:);

fid1 = fopen('kc.dat','w');
fid2 = fopen('pc.dat','w');

for i = 1:sk1+sk2+sk3
  fprintf(fid1,'%10.3f\n',kc(i,1));
end

for i = 1:sp1+sp2+sp3
  fprintf(fid2,'%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n',...
          pc(i,1),pc(i,2),pc(i,3),pc(i,4),pc(i,5),pc(i,6),pc(i,7));
end

fclose(fid1);
fclose(fid2);
