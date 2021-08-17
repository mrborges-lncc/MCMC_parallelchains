clear;
close all
loc=150;
jump=3;
N=00;
B=400;
A=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 5;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(1000, 1, Nchains);
Nfim = [2288 2379 2391 2499 2461 2369];
% Nfim = repmat(501, 1, Nchains);
Nfim = Nfim(Nch_ini+1:Nch_fim+1)-1;
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
nome = 'TwoPhase3D_onlyPerm_RW_RK';
nome = 'TwoPhase3D_RW_RK';
base_name = ['prod_D2_' nome];
hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
% hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMCrw_onlyPerm/';
% hom = '../';
% hom = '~/MCMCrw/';
homf= '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'twophaseflow/exp/prod/prod_referencia_0.dat']);
ref=dados;
home = [hom 'twoStage/select_prod/'];
%
file_name   = [home base_name '0_0.dat'];
data=load(file_name);
total=0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    pchains = load([home '../out/nchain_' nome n '.dat']);
    sz = size(pchains,1) + 1;
    total = total + sum(pchains(Nini(j)+1:Nfim(j)+1,2));
end
total= int64(total);
data = zeros(size(data,1),size(data,2),total);
%% MEDIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    name = [home '../out/nchain_' nome n '.dat']
    pchains = load(name);
    for i=Nini(j):Nfim(j)
        istr=num2str(i,5);
        file_name = [home base_name n '_' istr '.dat'];
        dat = load(file_name);
        m   = pchains(i+1,2);
        for nc = 1:m
            cont = cont + 1;
            data(:,:,cont) = dat;
        end
    end
end
fprintf('\n=========================================')
fprintf('\n=========================================')
fprintf('\nNumero total de dados: %d | %d',cont,total);
soma   = sum(data,2);
dmedio = mean(data,3);
erro   = std(data,0,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure1 = figure();

C=min(dados(:,1));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
D=max(dados(:,1))*1.01;

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados = dmedio;
errorbar(dados(1:jump:end,1),dados(1:jump:end,2),erro(1:jump:end,2),...
    'Parent',axes1,'Color',[0.85 0.33 0.10],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean 1','LineWidth',0.5);
errorbar(dados(1:jump:end,1),dados(1:jump:end,3),erro(1:jump:end,3),...
    'Parent',axes1,'Color',[0.07 0.62 1],'MarkerSize',4,'Marker','s',...
    'LineStyle','none','DisplayName','mean 2','LineWidth',0.5);
errorbar(dados(1:jump:end,1),dados(1:jump:end,4),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0.93 0.69 0.13],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean 3','LineWidth',0.5);
errorbar(dados(1:jump:end,1),dados(1:jump:end,5),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0 0 0],'MarkerSize',4,'Marker','s',...
    'LineStyle','none','DisplayName','mean 4','LineWidth',0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref(2:end,1),ref(2:end,2),'Parent',axes1,'Color',[0.85 0.33 0.10],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  1');
plot(ref(2:end,1),ref(2:end,3),'Parent',axes1,'Color',[0.07 0.62 1],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  2');
plot(ref(2:end,1),ref(2:end,4),'Parent',axes1,'Color',[0.93 0.69 0.13],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  3');
plot(ref(2:end,1),ref(2:end,5),'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
    'DisplayName','selection time');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(axes1,[0 D])
ylim(axes1,[A B])
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Oil rate ($m^3/day$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontName',...
    'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
    'latex','XMinorTick','on','YMinorTick','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast','FontSize',8);
set(legend1,'Box','off');

% Print
base=[homf base_name];
%print('-djpeg90',base)
print('-depsc','-r300',base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ERROM = mean(mean(erro(:,2:end)));
fprintf('\n=========================================')
fprintf('\n=========================================')
fprintf('\nDesvio padrao medio: %5.4f',ERROM)
fprintf('\n=========================================')
fprintf('\n=========================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Production
sz = size(data);
producao = zeros(sz(1),2,sz(3));
for k = 1:size(data,3)
    t = data(:,1,k);
    p = sum(data(:,2:end,k),2);
    pa= 0.0;
    prod = [0.0 0.0];
    for j = 2:size(t,1)
        dt = t(j) - t(j-1);
        pr = 0.5*(p(j) + p(j-1));
        pa = pa + pr * dt;
        prod = [prod; t(j) pa];
    end
    producao(:,:,k) = prod;
end
t = ref(:,1);
p = sum(ref(:,2:end),2);
prodref = [0.0 0.0];
pa= 0.0;
for j = 2:size(t,1)
    dt = t(j) - t(j-1);
    pr = 0.5*(p(j) + p(j-1));
    pa = pa + pr * dt;
    prodref = [prodref; t(j) pa];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
pause(1)
close all
figure2 = figure();

C=min(min(prodref(:,1,:)));
D=max(max(prodref(:,1,:)))*1.01;
A=min(min(prodref(:,2,:)));
B=max(max(prodref(:,2,:)))*1.01;

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure2,'FontSize',14,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');
% for k = 1:size(data,3)
%     h=plot(producao(1:end,1,k),producao(1:end,2,k),'Parent',axes1,...
%         'Color',[0.65 0.65 0.65],'MarkerSize',6,'LineWidth',0.5,...
%         'DisplayName','ref.  1');
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
plot(prodref(1:end,1),prodref(1:end,2),'Parent',axes1,...
    'Color',[1 0 0],'MarkerSize',6,'LineWidth',2,...
    'DisplayName','ref.')
prod = mean(producao,3);
plot(prod(1:end,1),prod(1:end,2),'Parent',axes1,...
    'Color',[0 0 0],'MarkerSize',6,'LineWidth',2,...
    'DisplayName','mean')
sprod = std(producao,0,3);
errorbar(prod(1:jump:end,1),prod(1:jump:end,2),sprod(1:jump:end,2),...
    'Parent',axes1,'Color',[0. 0. 0.],'MarkerSize',4,'Marker','none',...
    'LineStyle','-','DisplayName','error','LineWidth',0.5);

xlim(axes1,[C D]);
ylim(axes1,[A B]);
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Cumulative oil production ($m^3$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontName',...
    'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
    'latex','XMinorTick','on','YMinorTick','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',8);
set(legend1,'Box','off');

% Print
base=[homf 'total_' base_name];
%print('-djpeg90',base)
print('-depsc','-r300',base);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = min(size(ref,1),size(dados,1));
norma=norm(ref(1:N,2:end));
norma=norm(ref(1:N,2:end)-dados(1:N,2:end))/norma;
fprintf('ERRO RELATIVO = %e',norma);
fprintf('\n=========================================')
fprintf('\n=========================================')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
