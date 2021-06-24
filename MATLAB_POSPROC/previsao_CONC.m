clear;
close all
loc=100;
jump=4;
N=00;
B=450;
A=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 3;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(100, 1, Nchains);
Nfim = [130 156 182 105];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
base_name = 'prod_D2_KLfull_DE_RK'
dados=load('/home/mrborges/MCMCde/twophaseflow/exp/prod/prod_ref_0.dat');
ref=dados;
home = '/home/mrborges/MCMCde/twoStage/select_prod/'
%
data=zeros(size(ref,1),size(ref,2),sum(Nt));
k = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    for i=Nini(j):Nfim(j)
        k = k+1;
        istr=num2str(i,5);
        file_name   = [home base_name n '_' istr '.dat']
        data(:,:,k) = load(file_name);
    end
end
dmedio = mean(data,3);
erro   = std(data,0,3);
%
% Create figure
figure1 = figure()

C=min(dados(:,1));
D=max(dados(:,1))*1.01;

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');

dados = dmedio;
errorbar(dados(1:jump:end,1),dados(1:jump:end,2),erro(1:jump:end,2),...
    'Parent',axes1,'Color',[0.85 0.33 0.10],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean 1','LineWidth',1)
errorbar(dados(1:jump:end,1),dados(1:jump:end,3),erro(1:jump:end,3),...
    'Parent',axes1,'Color',[0.07 0.62 1],'MarkerSize',4,'Marker','s',...
    'LineStyle','none','DisplayName','mean 2','LineWidth',1)
errorbar(dados(1:jump:end,1),dados(1:jump:end,4),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0.93 0.69 0.13],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean 3','LineWidth',1)
errorbar(dados(1:jump:end,1),dados(1:jump:end,5),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0 0 0],'MarkerSize',4,'Marker','s',...
    'LineStyle','none','DisplayName','mean 4','LineWidth',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados=ref;
%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');
plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[0.85 0.33 0.10],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  1')
plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0.07 0.62 1],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  2')
plot(dados(:,1),dados(:,4),'Parent',axes1,'Color',[0.93 0.69 0.13],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  3')
plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
%     'DisplayName','time to select')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(axes1,[0 D])
ylim(axes1,[A B])
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Oil rate ($m^3/day$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast','FontSize',6);
set(legend1,'Box','off');

% Print
base=['./../../figuras/pres_' base_name];
%print('-djpeg90',base)
print('-depsc','-r100',base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma=norm(ref(:,2:end))
norma=norm(ref(:,2:end)-dados(:,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
