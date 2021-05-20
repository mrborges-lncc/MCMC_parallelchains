clear;
dados=load('../prod/prodF_ref_0.dat');
%dados=load('../conc/conc_amostra_0.dat');
%dados=load('../../MCMC/select_prod/prodF_3.dat');
%dados=load('../../MCMC/reject_prod/prodF_30.dat');
%
% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1);
%box(axes1,'on');

plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],...
    'Marker','o','LineStyle','none','DisplayName','ref1')
hold on
plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],...
    'Marker','o','LineStyle','none','DisplayName','ref2')
% plot(dados(:,1),dados(:,4),'Color',[1 0 0],'Marker','^','LineStyle','none')
% plot(dados(:,1),dados(:,5),'Color',[0 0 1],'Marker','+','LineStyle','none')
% plot(dados(:,1),dados(:,6),'Color',[0 1 0],'Marker','v','LineStyle','none')
% plot(dados(:,1),dados(:,7),'Color',[0 0 0],'Marker','s','LineStyle','none')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados=load('../../MCMC/select_prod/prodF_18.dat');
%dados=load('../../MCMC/reject_prod/prodF_0.dat');

plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],'LineWidth',2,'DisplayName','amostra1')
plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],'LineWidth',2,'DisplayName','amostra2')
% plot(dados(:,1),dados(:,4),'Color',[1 0 0],'LineWidth',2)
% plot(dados(:,1),dados(:,5),'Color',[0 0 1],'LineWidth',2)
% plot(dados(:,1),dados(:,6),'Color',[0 1 0],'LineWidth',2)
% plot(dados(:,1),dados(:,7),'Color',[0 0 0],'LineWidth',2)

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','SouthWest');
clear
