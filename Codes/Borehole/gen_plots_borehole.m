close all;
clear all;
addpath(genpath('./line_fewer_markers_v4/'));

err_samples();
%ub_conv();

function ubc = ub_conv()
data_ub = load('nsams_ub.txt');
%data_sti = load('STi_conv.mat');
col = ['c','r','b','m','g','k'];
ls = ['-','-.','--','-s','-d',':','-*'];

figure
hold on;
%for i = 1:7
%plot(data_ub(:,1),data_ub(:,i+1),'LineWidth',2,'color',col(i));
%plot(data_ub(:,1),data_ub(:,2),'LineWidth',1.5,'color',[0 0 0]);
%plot(data_ub(:,1),data_ub(:,3),':','LineWidth',1.5,'color',[0 0 0]);
%plot(data_ub(:,1),data_ub(:,4),'-.','LineWidth',1.5,'color',[0 0 0]);
%plot(data_ub(:,1),data_ub(:,5),':','LineWidth',1.5,'color',[0 0 0]);
%plot(data_ub(:,1),data_ub(:,6),'-.','LineWidth',1.5,'color',[0 0 0]);
%plot(data_ub(:,1),data_ub(:,7),'-.','LineWidth',1.5,'color',[0 0 0]);
line_fewer_markers(data_ub(:,1),data_ub(:,2),20,':ko','MarkerSize',6,'MarkerFaceColor','k');
plot(data_ub(:,1),data_ub(:,3),'--k');
line_fewer_markers(data_ub(:,1),data_ub(:,4),20,':ks','MarkerSize',6,'MarkerFaceColor','k');
plot(data_ub(:,1),data_ub(:,5),'--k');
line_fewer_markers(data_ub(:,1),data_ub(:,6),20,':ks','MarkerSize',6,'MarkerFaceColor','k');
%line_fewer_markers(data_ub(:,1),data_ub(:,7),10,':kd','MarkerSize',4,'MarkerFaceColor','k');
plot(data_ub(:,1),data_ub(:,7),'-.k');
line_fewer_markers(data_ub(:,1),data_ub(:,8),20,':k*','MarkerSize',6);
%line_fewer_markers(t, sin(3*t).*cos(t/2), 10, 'p-');
%end

hold off;
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{\hat{\mathcal{C}_i\mu_i}}$$','interpreter','latex');
ylim([-0.01,1.1.*max(max(data_ub(:,2:end)))]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{r_w}~[0.6941]$','$\mathrm{T_u}~[0.0000]$',...
              '$\mathrm{H_u}~[0.1061]$','$\mathrm{T_l}~[0.0000]$',...
              '$\mathrm{H_l}~[0.1061]$','$\mathrm{L}~[0.1028]$',...
              '$\mathrm{K_w}~[0.0251]$'},'location','East');
set(leg,'Interpreter','latex');
box on;
print -depsc ub_conv_borehole.eps
end

function err = err_samples()
%ns = [10;20;30;40;50;60;70;80];
ns = [10;20;30;40;50];
%err7D = [1.836;2.958e-1;1.864e-2;1.632e-2;1.934e-2;2.374e-3;2.319e-3;1.184e-3];
%err5D = [4.018e-1;3.136e-1;1.629e-2;2.955e-2;2.345e-3;2.770e-3;4.879e-4;4.331e-4];
%err4D = [1.135e-1;1.603e-1;7.151e-3;1.571e-3;1.519e-4;5.217e-5;3.378e-4;4.286e-6];
err7D = [1.836;2.958e-1;1.864e-2;1.632e-2;1.934e-2];
err5D = [4.018e-1;3.136e-1;1.629e-2;2.955e-2;2.345e-3];
err4D = [1.135e-1;1.603e-1;7.151e-3;1.571e-3;1.519e-4];
err = zeros(size(ns,1),3);
err(:,1) = log10(err7D(:));
err(:,2) = log10(err5D(:));
err(:,3) = log10(err4D(:));

figure
hold on;
plot(ns,err(:,1),':ko','MarkerFaceColor','k');
plot(ns,err(:,2),':kd','MarkerFaceColor','k');
plot(ns,err(:,3),':ks','MarkerFaceColor','k');
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{log_{10}(\epsilon_{LOO})}$$','interpreter','latex');
xlim([8,52]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{7D PCE}$','$\mathrm{5D PCE}$','$\mathrm{4D PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');
box on;
print -depsc err_samples_borehole.eps
end



