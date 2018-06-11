close all
clear all

%plot_ub();

%function pub = plot_ub()
ub5 = load('new_lean_mix/ub5.mat');
ub10 = load('new_lean_mix/ub10.mat');
ub15 = load('new_lean_mix/ub15.mat');
ub20 = load('new_lean_mix/ub20.mat');
mu5 = load('new_lean_mix/mu5.mat');
mu10 = load('new_lean_mix/mu10.mat');
mu15 = load('new_lean_mix/mu15.mat');
mu20 = load('new_lean_mix/mu20.mat');
xp = linspace(1,19,19);

figure;
hold on;
h = plot(xp,ub5.ub,'--o','MarkerFaceColor','r');
plot(xp,ub10.ub,'--o','MarkerFaceColor','b');
plot(xp,ub15.ub,'--o','MarkerFaceColor','g');
plot(xp,ub20.ub,'--o','MarkerFaceColor','k');
xticklabels({'$\mathrm{A_1}$','$\mathrm{A_2}$','$\mathrm{A_3}$','$\mathrm{A_4}$',...
             '$\mathrm{A_5}$','$\mathrm{A_6}$','$\mathrm{A_7}$','$\mathrm{A_8}$',...
             '$\mathrm{A_9}$','$\mathrm{A_{10}}$','$\mathrm{A_{11}}$','$\mathrm{A_{12}}$',...
             '$\mathrm{A_{13}}$','$\mathrm{A_{14}}$','$\mathrm{A_{15}}$','$\mathrm{A_{16}}$',...
             '$\mathrm{A_{17}}$','$\mathrm{A_{18}}$','$\mathrm{A_{19}}$'});
set(gca,'xtick',1:19,'fontsize',10,'TickLabelInterpreter','latex');
ax = ancestor(h, 'axes');
yrule = ax.YAxis;
yrule.FontSize = 16;
xlabel('$$\mathrm{Parameter}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\widehat{\mathcal{C}_i\mu_i}}$$','interpreter','latex','fontsize',20);
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{N~=~5}$','$\mathrm{N~=~10}$','$\mathrm{N~=~15}$',...
              '$\mathrm{N~=~20}$','location','best'});
set(leg,'Interpreter','latex','fontsize',14);
box on;
grid on;
print -depsc ub_conv_kinetics_lean.eps

figure;
mu_max = zeros(3,1);
mu_max(1,1) = max(abs(mu5.mu - mu10.mu));
mu_max(2,1) = max(abs(mu10.mu - mu15.mu));
mu_max(3,1) = max(abs(mu15.mu - mu20.mu));
%yl = 10e-17./3.16;
%yu = 10e-17.*3.16;
semilogy(1:3,mu_max,'--o','MarkerFaceColor','k');
xlabel('$$\mathrm{Iteration}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\Delta\mu_s}$$','interpreter','latex','fontsize',20);
xlim([0,4]);
%ylim([yl,yu]);
set(gcf,'color',[1,1,1]);
set(gca,'xtick',1:3,'fontsize',18,'TickLabelInterpreter','latex');
box on;
grid on;
print -depsc mu_lean.eps
%end



