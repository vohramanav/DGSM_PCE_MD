close all;
clear all;
rng(50);

dim = 2;
L = zeros(1,dim); U = zeros(1,dim);

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1) = 0.9.*nom(1); U(1) = 1.1.*nom(1);
L(2) = 0.9.*nom(9); U(2) = 1.1.*nom(9);

nsams = 100; pts = zeros(nsams,19);

x = linspace(L(1),U(1),10);
y = linspace(L(2),U(2),10);

[X,Y] = meshgrid(x,y);

f1 = fopen('pts_smooth_rich.txt','w');
fmt = '%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n';

for j = 2:8
  pts(:,j) = nom(j);
end

for j = 10:19
  pts(:,j) = nom(j);
end

c = 1;
for j = 1:10
  for i = 1:10
   pts(c,1) = X(j,i);
   pts(c,9) = Y(j,i);
   fprintf(f1,fmt,pts(c,1),pts(c,2),pts(c,3),pts(c,4),pts(c,5),pts(c,6),pts(c,7),pts(c,8),pts(c,9),pts(c,10),pts(c,11),pts(c,12),pts(c,13),pts(c,14),pts(c,15),pts(c,16),pts(c,17),pts(c,18),pts(c,19));
   c = c + 1;    
  end
end

fclose(f1);

ig = load('record_id_smooth_rich.txt');
ig_surf = reshape(ig,size(X));
save('X.mat','X');
save('Y.mat','Y');
save('ig_surf.mat','ig_surf');

figure;
pcolor(X,Y,ig_surf);
shading interp;
%surf(X,Y,ig_surf);
c = colorbar();
xlabel('$$\mathrm{A_1}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{A_9}$$','interpreter','latex','fontsize',20);
set(gca,'TickLabelInterpreter','latex','fontsize',18);
set(c,'TickLabelInterpreter','latex')
set(gcf,'color',[1,1,1]);
print -dpng ig_rich.png
