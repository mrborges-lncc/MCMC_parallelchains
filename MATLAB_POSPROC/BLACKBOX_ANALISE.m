clear;
barrel=1;%.158987295; %m^3
base_name = 'prod_D1_blackbox_RK7_950';
base_name1 = 'prod_D2_blackbox_RK7_950';
base_name2 = 'prod_D3_blackbox_RK7_950';
%dados=load('../prod/prodF_ref_0.dat');
%dados=load('../simuladorBTMM/exp/prod/prodF_ref_0.dat');
dados=load('../blackbox/exp/saidaref1.dat');
dados1=load('../blackbox/exp/saidaref2.dat');
dados2=load('../blackbox/exp/saidaref3.dat');
%dados=load('../simul_comp/exp/prod/prodF_amostra_0.dat');
%dados=load('../twoStage/simuladorBTMM/exp01/prod/prodF_ref_0.dat');
%dados=load('../../MCMC/select_prod/prodF_3.dat');
%dados=load('../../MCMC/reject_prod/prodF_30.dat');
dados(:,2)=dados(:,2)/barrel;
dados1(:,2)=dados1(:,2)/barrel;
dados2(:,2)=dados2(:,2)/barrel;
ref =dados;
ref1=dados1;
ref2=dados2;
%
% Create figure
% Create figure
figure1 = figure()
% B=15;
% A=0;
A=min(min(ref2(:,2)),min(min(ref1(:,2)),min(ref(:,2))))
if(A<0)
    A=A*1.2;
else
    A=A*0.5;
end
B=max(max(dados2(:,2)),max(max(dados1(:,2)),max(dados(:,2))))*1.1;
C=min(min(dados2(:,1)),min(max(dados1(:,1)),min(dados(:,1))))*1.1;
D=max(max(dados2(:,1)),max(max(dados1(:,1)),max(dados(:,1))))*1.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name = ['../twoStage/select_prod/' base_name '.dat']
dados=load(file_name);
dados(:,2)=dados(:,2)/barrel;
file_name1 = ['../twoStage/select_prod/' base_name1 '.dat']
dados1=load(file_name1);
dados1(:,2)=dados1(:,2)/barrel;
file_name2 = ['../twoStage/select_prod/' base_name2 '.dat']
dados2=load(file_name2);
dados2(:,2)=dados2(:,2)/barrel;

dasp=[1 1.*(B-A)/(D-C) 1];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');


plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],...
    'MarkerSize',8,'Marker','o','LineStyle','none','DisplayName','sample 1')
plot(dados1(:,1),dados1(:,2),'Parent',axes1,'Color',[0 0 1],...
    'MarkerSize',8,'Marker','s','LineStyle','none','DisplayName','sample 2')
plot(dados2(:,1),dados2(:,2),'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',8,'Marker','^','LineStyle','none','DisplayName','sample 3')
% plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',5,'Marker','+','LineStyle','none','DisplayName','sample 4')
% plot(dados(:,1),dados(:,6),'Parent',axes1,'Color',[0 1 0],...
%     'MarkerSize',5,'Marker','v','LineStyle','none','DisplayName','sample 5')
% plot(dados(:,1),dados(:,7),'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',5,'Marker','s','LineStyle','none','DisplayName','sample 6')
%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');
%dados=load('../twoStage/simuladorBTMM/exp01/prod/prodF_ref_0.dat');

plot(ref(:,1),ref(:,2),'Parent',axes1,'Color',[1 0 0],...
    'LineWidth',2,'DisplayName','ref. 1')
plot(ref1(:,1),ref1(:,2),'Parent',axes1,'Color',[0 0 1],...
    'LineWidth',2,'DisplayName','ref. 2')
plot(ref2(:,1),ref2(:,2),'Parent',axes1,'Color',[0 0 0],...
    'LineWidth',2,'DisplayName','ref. 3')
% plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
%     'LineWidth',2,'DisplayName','ref. 4')
% plot(dados(:,1),dados(:,6),'Parent',axes1,'Color',[0 1 0],...
%     'LineWidth',2,'DisplayName','ref. 5')
% plot(dados(:,1),dados(:,7),'Parent',axes1,'Color',[0 0 0],...
%     'LineWidth',2,'DisplayName','ref. 6')

% Create xlabel
xlim(axes1,[C D])
ylim(axes1,[A B])
% Create xlabel
xlabel('t','FontSize',20,'FontName','Times New Roman','FontAngle','italic');

% Create ylabel
ylabel('F(STB/day)','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',12);
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
