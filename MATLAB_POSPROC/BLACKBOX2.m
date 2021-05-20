clear;
barrel=1;%.158987295; %m^3
base_name = 'prod_D3_blackbox_RK7_750';
%dados=load('../prod/prodF_ref_0.dat');
%dados=load('../simuladorBTMM/exp/prod/prodF_ref_0.dat');
dados=load('../blackbox/exp/saidaref3.dat');
%dados=load('../simul_comp/exp/prod/prodF_amostra_0.dat');
%dados=load('../twoStage/simuladorBTMM/exp01/prod/prodF_ref_0.dat');
%dados=load('../../MCMC/select_prod/prodF_3.dat');
%dados=load('../../MCMC/reject_prod/prodF_30.dat');
dados(:,2)=dados(:,2)/barrel;
ref=dados;
%
% Create figure
% Create figure
figure1 = figure()
B=15;
A=0;
A=min(dados(:,2));
if(A<0)
    A=A*1.2;
else
    A=A*0.9;
end
B=max(dados(:,2))*1.1;
C=min(dados(:,1));
D=max(dados(:,1))*1.;

dasp=[1 1.*(B-A)/(D-C) 1];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');


plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],...
    'MarkerSize',5,'Marker','o','LineStyle','none','DisplayName','ref. 1')
% plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],...
%     'MarkerSize',5,'Marker','o','LineStyle','none','DisplayName','ref. 2')
% plot(dados(:,1),dados(:,4),'Parent',axes1,'Color',[1 0 0],...
%     'MarkerSize',5,'Marker','^','LineStyle','none','DisplayName','ref. 3')
% plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',5,'Marker','+','LineStyle','none','DisplayName','ref. 4')
% plot(dados(:,1),dados(:,6),'Parent',axes1,'Color',[0 1 0],...
%     'MarkerSize',5,'Marker','v','LineStyle','none','DisplayName','ref. 5')
% plot(dados(:,1),dados(:,7),'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',5,'Marker','s','LineStyle','none','DisplayName','ref. 6')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name = ['../twoStage/select_prod/' base_name '.dat']
dados=load(file_name);
dados(:,2)=dados(:,2)/barrel;
%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');
%dados=load('../twoStage/simuladorBTMM/exp01/prod/prodF_ref_0.dat');

plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],...
    'LineWidth',2,'DisplayName','simul. 1')
% plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],...
%     'LineWidth',2,'DisplayName','simul. 2')
% plot(dados(:,1),dados(:,4),'Parent',axes1,'Color',[1 0 0],...
%     'LineWidth',2,'DisplayName','simul. 3')
% plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
%     'LineWidth',2,'DisplayName','simul. 4')
% plot(dados(:,1),dados(:,6),'Parent',axes1,'Color',[0 1 0],...
%     'LineWidth',2,'DisplayName','simul. 5')
% plot(dados(:,1),dados(:,7),'Parent',axes1,'Color',[0 0 0],...
%     'LineWidth',2,'DisplayName','simul. 6')

% Create xlabel
xlim(axes1,[C D])
ylim(axes1,[A B])
% Create xlabel
xlabel('t','FontSize',20,'FontName','Times New Roman','FontAngle','italic');

% Create ylabel
ylabel('F(STB/day)','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast','FontSize',12);
set(legend1,'Box','off');

% Print
base=['../figuras/conc_' base_name];
print('-djpeg90',base)
print('-depsc','-r300',base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma=norm(ref(:,2:end))
norma=norm(ref(:,2:end)-dados(:,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
