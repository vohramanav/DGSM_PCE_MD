close all
clear all
addpath(genpath('/home/manavvohra/jw_tools/'));
warning('off','all');

%sq = load('samples_qoi.txt');
q1 = 0.0;
q2 = 0.1;
g1 = 1.08;
g2 = 1.32;

ndim=2;
nord=4;
Nq1d=nord+1;

multiIndex = mult_ind_iso(nord,ndim);
nPCTerms = size(multiIndex,1);
[X,w] = gen_full_quad(repmat(Nq1d,1,ndim),'GL');

nisp = nisp_gen_xw(multiIndex,X,w);

Nq = size(w,1);
x_nodes = X(1:Nq1d,1);

% %______________ REALIZATIONS ______________%

Xr = zeros(size(X,1),2);
Xr(:,1) = q1 + ((q2-q1)./2).*(X(:,1) + 1);
Xr(:,2) = g1 + ((g2-g1)./2).*(X(:,2) + 1);

%F = sq(:,3);
%Ck = nisp*F;

% %______________ ERROR ___________________%

% Project real Sobol samples in the interval [-1,1] 

%Xsobol = 2.0.*(sobol(:,1)-l1)./(l2-l1) - 1.0;
%Ysobol = 2.0.*(sobol(:,2)-g1)./(g2-g1) - 1.0;
%pce_sobol = reconstructPC([Xsobol Ysobol],Ck,multiIndex);
%NR = sum(((sobol(:,3)-pce_sobol).^2)).^(0.5);
%DR = sum((sobol(:,3).^2)).^(0.5);
%E = NR./DR;

% %______________ GSA ______________________%

%TS = zeros(2,1);
%FS = zeros(2,1);
%[TS(:),~,FS(:)]=sensitivity_calc(Ck,multiIndex);
%
%for i = 1:2
%    ST(i,1) = FS(i);
%    ST(i,2) = TS(i);
%end



% %________________ PLOTS __________________%

% Plotting error surface
% figure
% [XX,YY] = meshgrid(linspace(-1,1,40));
% data_pce = reconstructPC([XX(:) YY(:)],Ck,multiIndex);
% data_pce_surf = reshape(data_pce,size(XX));
% XXr = l1 + ((l2-l1)./2).*(XX + 1);
% YYr = g1 + ((g2-g1)./2).*(YY + 1);
% pcolor(XXr,YYr,data_pce_surf);
% c = colorbar();
% hold on
% p2=scatter(Xr(:,1),Xr(:,2),'filled','MarkerFaceColor','k');
% %p3=scatter(sobol(:,1),sobol(:,2),'filled','d','MarkerFaceColor','r');
% xlabel('$$\mathrm{L~(\AA)}$$','interpreter','latex','fontsize',18)
% ylabel('$$\mathrm{dT/dx~(\frac{K}{\AA})}$$','interpreter','latex','fontsize',18)
% shading interp;
% leg = legend([p2],{'$\mathrm{G-L~Nodes}$'},'location','northoutside');
% %leg = legend([p3],{'$\mathrm{Sobol~Samples}$'},'location','northoutside');
% set(leg,'Interpreter','latex')
% set(gca,'fontsize',14);
% set(gca,'TickLabelInterpreter','latex')
% set(c,'TickLabelInterpreter','latex')
% set(gcf,'color',[1,1,1]);
% print -dpng err2D_300.png
%
%% Plotting the spectrum
%num = linspace(1,15,15);
%hline = linspace(0,0,15);
%sz = 50.*(log10(abs(Ck))+1.1*max(abs(log10(Ck))));
%c = zeros(15,1);
%figure
%p1 = scatter(num(1),log10(abs(Ck(1))),sz(1),'g','filled');
%hold on
%p2 = scatter(num(2),log10(abs(Ck(2))),sz(2),'b','filled');
%p3 = scatter(num(3:end),log10(abs(Ck(3:end))),sz(3:end),'m','filled');
%p4 = plot(num,hline,'--','color','k');
%leg = legend([p1 p2 p3],{'$\mathrm{Mean}$','$\mathrm{p~=~1}$','$\mathrm{p~>~2}$'});
%xlim([0,16]);
%xlabel('$$\mathrm{k}$$','interpreter','latex','fontsize',18);
%ylabel('$$\mathrm{log_{10}(|c_{k}|)}$$','interpreter','latex','fontsize',18);
%set(leg,'Interpreter','latex');
%set(gca,'fontsize',14);
%set(gca,'TickLabelInterpreter','latex')
%set(gcf,'color',[1,1,1]);
%box on
%print -dpng PCspectrum_300.png

% % Sensitivity Plot

%fig = figure;
%%bar(Yfirst,'stacked','Linestyle','none');
%bar(ST);
%ylabel('Sensitivity','fontsize',16);
%set(gca,'XTickLabel',{'L','dT/dx'});
%%title('T = 900 K','fontsize',16);
%set(gca,'fontsize',16);
%set(gcf,'color',[1 1 1]);
%colormap linspecer;





