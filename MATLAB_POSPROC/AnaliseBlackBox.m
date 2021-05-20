close all;
clear;
analysis = 1;
convanalysis = 1;
tempo = cputime;
d  = 4; % dimension of multivariate Gausesian
Nc0= 0; % Numero da cadeia inicial
Nc = 15; % Numero da cadeia final
NT = Nc-Nc0+1; % Numero total de cadeias
NCAD = [Nc0:1:Nc];
MM   = 500;
transp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
tm = [];
rep=[];
home = '../coarrayMCMC/';
% home = '../twoStage/';
base_name = 'blackboxRW_RK';
%base_name = 'blackboxDE_RK';
% base_name = 'blackboxAM_RK';
% base_name = 'blackboxDREAM_1X1_RK';
nome_extra = '';
for i=Nc0:Nc
    file_name = ...
        [home 'error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
    Ncurv = size(dados,2)-2;
    file_name =...
        [home 'out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    dados=[dados; sz 1];
    rep = [rep;dados];
end
tm=tm-1;
Nfim = tm-1;
Nini =0*Nfim+MM;
N=max(tm)+1;
clear dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DADOS DE REFERENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
refile = '../blackbox/exp/saidaref';
szref  = 1;
ref    = [];
for i=1:Ncurv
    nomef = [refile num2str(i,3) '.dat']
    dataref = load(nomef);
    sz = size(dataref,1);
    szref = max(szref,sz);
    ref = [ref dataref];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variavel  = [home 'select_thetas/theta_v'];
variavel1 = [home 'out/nchain_'];
variavel2 = [home 'select_prod/prod_D'];
nchain=zeros(N,NT);
nk=0;
SOMA = [];
mu1 = zeros(1,d);
mu2 = zeros(1,d);
mu3 = zeros(1,d);
m2u1 = zeros(1,d);
m2u2 = zeros(1,d);
m2u3 = zeros(1,d);
media= zeros(3,d);
scur1= zeros(sz,2);
scur2= zeros(sz,2);
scur3= zeros(sz,2);
dcur1= zeros(sz,2);
dcur2= zeros(sz,2);
dcur3= zeros(sz,2);
nk = 0;
for i=1:NT
    nfile =...
        [variavel1 base_name num2str(NCAD(i),3) '.dat'];
    rep = load(nfile);
%     nchain = 1+0*rep(:,2);
    nchain = rep(:,2);
    for k=Nini(i):Nfim(i)
        nfile1 =...
            [variavel '1_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        nfile2 =...
            [variavel '2_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        nfile3 =...
            [variavel '3_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        nk=nk+nchain(k);
        x1 = load(nfile1)';
        x2 = load(nfile2)';
        x3 = load(nfile3)';
        nfile1 =...
            [variavel2 '1_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        nfile2 =...
            [variavel2 '2_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        nfile3 =...
            [variavel2 '3_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        curva1 = load(nfile1);
        curva2 = load(nfile2);
        curva3 = load(nfile3);
        scur1  = scur1+nchain(k)*curva1;
        scur2  = scur2+nchain(k)*curva2;
        scur3  = scur3+nchain(k)*curva3;
        dcur1  = dcur1+nchain(k)*curva1.^2;
        dcur2  = dcur2+nchain(k)*curva2.^2;
        dcur3  = dcur3+nchain(k)*curva3.^2;
        mu1= mu1 + nchain(k)*x1;
        mu2= mu2 + nchain(k)*x2;
        mu3= mu3 + nchain(k)*x3;
        m2u1= m2u1 + nchain(k)*x1.^2;
        m2u2= m2u2 + nchain(k)*x2.^2;
        m2u3= m2u3 + nchain(k)*x3.^2;
    end
end
mu1 = mu1/nk;
mu2 = mu2/nk;
mu3 = mu3/nk;
m2u1 = m2u1/nk;
m2u2 = m2u2/nk;
m2u3 = m2u3/nk;
scur1= scur1/nk;
scur2= scur2/nk;
scur3= scur3/nk;
dcur1= dcur1/nk;
dcur2= dcur2/nk;
dcur3= dcur3/nk;
sd1  = sqrt(dcur1(:,2)-scur1(:,2).^2);
sd2  = sqrt(dcur2(:,2)-scur2(:,2).^2);
sd3  = sqrt(dcur3(:,2)-scur3(:,2).^2);
fprintf('\n Desvios médios: sigma_1 = %f; sigma_2 = %f; sigma_3 = %f \n',...
    mean(sd1),mean(sd2),mean(sd3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info=[];
for i=Nc0:Nc
    file_name = [home 'error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
    data =[data; dados];
    file_name = [home 'out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    dados=[dados; sz 1];
    rep = [rep;dados];
    file_name = [home 'in/init_stat_' base_name num2str(i,'%1.1d') '.in']
    fileID = fopen(file_name);
    [A] = fscanf(fileID,'%d %d %d %d');
    info= [info; A(1) A(3)];
    fclose(fileID);
end
taxa=100.0*double(info(:,2))./double(info(:,1));
AR = mean(taxa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = max(max(max(ref(:,2)),max(ref(:,4))),max(ref(:,6)))*1.2;
A = min(min(min(ref(:,2)),min(ref(:,4))),min(ref(:,6)))*1.2;
A = -100;
B = 100;
xdata = scur1(:,1);
% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',18,...
    'FontName','Times New Roman',...
    'DataAspectRatio',[1 1*(B-A)/(max(xdata)*1.-min(xdata)*1.) 1]);
box(axes1,'on');
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref(:,1),ref(:,2),'r-','LineWidth',3);
hold on
plot(ref(:,3),ref(:,4),'b-','LineWidth',3);
plot(ref(:,5),ref(:,6),'k-','LineWidth',3);
errorbar(scur1(:,1),scur1(:,2),sd1,'Marker','^','MarkerSize',4,...
    'LineWidth',2,'Color',[1 0 0],'LineStyle','none')
errorbar(scur2(:,1),scur2(:,2),sd2,'Marker','s','MarkerSize',4,...
    'LineWidth',2,'Color',[0 0 1],'LineStyle','none')
errorbar(scur3(:,1),scur3(:,2),sd3,'Marker','o','MarkerSize',4,...
    'LineWidth',2,'Color',[0 0 0],'LineStyle','none')
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[min(xdata)*1.1 max(xdata)*1.1]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[A B]);
box(axes1,'on');
hold(axes1,'all');
% Create xlabel
xlabel('x','FontSize',18,'FontName','Times New Roman','FontAngle','italic');
% Create ylabel
ylabel('f(x)','FontSize',18,'FontName','Times New Roman','FontAngle','italic');
% Create textbox
local = [0.22 0.15 0.18 0.067];
local = [0.65 0.8 0.18 0.067];
local = [0.6 0.8 0.26 0.067];
annotation(figure1,'textbox',...
    local,...
    'String',{['AR = ' num2str(AR,'%4.1f') ' %']},...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Times',...
    'FitBoxToText','off',...
    'LineStyle','none');
nome = '../figuras/Blackbox_MCMCpost_';
name = [nome base_name '_' nome_extra]
print('-depsc','-r300',name)
s1 = (mean(sd1))^2;
s2 = (mean(sd2))^2;
s3 = (mean(sd3))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(convanalysis == 1)
    N = info(1,1);
    x = zeros(1,d);
    for i=Nc0:Nc
        file_name = ...
            [home 'error/erros_' base_name num2str(i,'%1.1d') '.dat'];
        dados=load(file_name);
        sz = size(dados,1);
        tm = [tm;sz];
        Ncurv = size(dados,2)-2;
        file_name =...
            [home 'out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
        dados=load(file_name);
        dados=[dados; sz 1];
        rep  = [rep;dados];
    end
    for nd=1:Ncurv
        X = zeros(N,d,NT);
        for i=1:NT
            nel=0;
            nfile =...
                [variavel1 base_name num2str(NCAD(i),3) '.dat']
            rep = load(nfile);
            nchain = rep(:,2);
            nchain = [nchain; N-sum(nchain)];
            for k=0*Nini(i):Nfim(i)+1
                nfile =...
                    [variavel num2str(nd,3) '_' base_name...
                    num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
                x = load(nfile)';
                for nrep=1:nchain(k+1)
                    nel = nel+1;
                    X(nel,:,i) = x;
                end
            end
        end
%% convergencia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cv1 = [];
        xcv = [];
        V1  = [];
        W1  = [];
        b = N/20;
        for i=b:b:N
            ini = i/2;
            fim = i;
            xcv = [xcv; fim];
%             [cv,v,w] = convergeR(X(ini:fim,:,:));
            [cv,v,w] = convergeRMatrix(X(ini:fim,:,:));
            cv1 = [cv1; cv];
            W1  = [W1; w];
            V1  = [V1; v];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rfigure(xcv,cv1,transp)
        nome = '../figuras/Blackbox_MCMC_R_';
        name = [nome base_name '_' nome_extra 'v' num2str(nd,3)]
        print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rfigure2(xcv,(W1),(V1),transp)
        nome = '../figuras/Blackbox_MCMC_WV_';
        name = [nome base_name '_' nome_extra 'v' num2str(nd,3)]
        print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
fprintf('\n Desvios médios: sigma_1 = %f; sigma_2 = %f; sigma_3 = %f \n',...
    mean(sd1),mean(sd2),mean(sd3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
