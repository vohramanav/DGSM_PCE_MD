close all
clear all

ub_conv();
%err_samples();

function ubc = ub_conv()
data_ub = load('nsams_ub.txt');
data_sti = load('STi_conv.mat');
col = ['c','r','b','m','g','k'];

figure
hold on;
for i = 1:size(col,2)
plot(data_ub(:,1),data_ub(:,i+1),'LineWidth',2,'color',col(i));
end

hold off;
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{\hat{ub_i}}$$','interpreter','latex');
ylim([-0.01,1.1.*max(max(data_ub(:,2:end)))]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{m}~[0.0051]$','$\mathrm{k_1}~[0.0271]$',...
              '$\mathrm{k_2}~[0.0003]$','$\mathrm{r}~[0.2596]$',...
              '$\mathrm{F}~[0.3952]$','$\mathrm{t_1}~[0.3279]$'},...
              'location','northeast');
set(leg,'Interpreter','latex');
box on;
print -depsc ub_conv_oscillator.eps
end

function err = err_samples()
ns = [10;20;30;40;50];
err6D = [1.588;3.197e-2;1.275e-3;5.584e-4;1.573e-3];
err4D = [1.616e-1;4.063e-3;9.928e-4;2.943e-5;1.747e-5];
err3D = [8.999e-4;3.688e-5;1.035e-8;4.155e-12;7.728e-12];
err = zeros(size(ns,1),3);
err(:,1) = log10(err6D(:));
err(:,2) = log10(err4D(:));
err(:,3) = log10(err3D(:));
xp = 8:1:52;
yp = log10(1e-6);

figure
hold on;
plot(ns,err(:,1),'-o','color','k','linewidth',2);
plot(ns,err(:,2),'-d','color','b','linewidth',2);
plot(ns,err(:,3),'-s','color','r','linewidth',2);
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{log_{10}(\epsilon_{LOO})}$$','interpreter','latex');
xlim([8,52]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{6D PCE}$','$\mathrm{4D PCE}$','$\mathrm{3D PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');
box on;
grid on;
print -depsc err_samples_oscillator.eps
end
