clear;
close all
loc=100;
jump=1;
jumpstat=1;
P=0.90; % probability
intconf= 10; % if == 1 intervalos de confianca cc. +/- desvio*t
printa = 1; % if == 1 imprime os intervalos
printa2= 10; % if == 1 imprime as amostras
N=21;
Nini=50;
Nfim=99;
Nt=(Nfim-Nini)+1;
base_name = 'prodF_MC_';
base_name = 'prodF_MCMC_';
base_name = 'blackbox_RK7_';
ref=load('../blackbox/exp/saidaref1.dat');
variavel = '../twoStage/select_prod/prod_D1_';
X=[loc loc];
x1=[-1 1];
%y1=[-3 20];
%y=[min(ref(:,2)) max(ref(:,2))];
yl=[(min(min(ref(:,2:end)))-max(max(ref(:,2:end)))*0.25) max(max(ref(:,2:end)))*1.25];
%yl=[0.1 max(max(ref(:,2:end)))];
xl=[min(min(ref(:,1))) max(max(ref(:,1)))];
vy=abs(min(min(ref(:,2:end))) - max(max(ref(:,2:end))));
vx=abs(min(ref(:,1)) - max(ref(:,1)));
vy=yl(2)-yl(1);
%yl=y1
%
% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
    'DataAspectRatio',[1 1*vy/vx 1]);
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min2=1e20;
for i=Nini:Nfim
    istr=num2str(i,5);
    file_name = [variavel base_name istr '.dat']
    data=load(file_name);
    dados=data(1:N,:);
    aux=min(min(dados(:,2:end)));
    min2=min(min2,aux);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chegada=zeros(Nt,1);
pt=0;
DATA=zeros(N,size(ref,2)-1,Nt);
m=0;
for i=Nini:Nfim
    m=1+m;
    istr=num2str(i,5);
    file_name = [variavel base_name istr '.dat'];
    data=load(file_name);
    dados=data(1:N,:);
    DATA(:,:,m)=dados(:,2:end);
    for j=2:size(dados,2)
        cor = [0.8 0.8 0.8];
        if(j==3)
            cor = [0.5 0.5 0.5];
        end
        if(printa2==1)
            plot(dados(:,1),dados(:,j),'Parent',axes1,...
                'MarkerEdgeColor',[0.8314 0.8157 0.7843],'MarkerSize',4,...
                'Marker','none','LineWidth',1,'Color',cor);
        end
    end
end
Dmedia=mean(DATA,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTERVALO DE CONFIANCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sdesv=zeros(N,size(dados,2)-1);
for j=1:size(dados,2)-1
    for i=1:N
        sdesv(i,j)=sqrt(var(DATA(i,j,:)));
    end
end
if(intconf==1)
    range=tinv((1-(1-P)/2),Nt-1)*sdesv/sqrt(Nt);
    intervS=Dmedia+range;
    intervI=Dmedia-range;
else
    range=tinv((1-(1-P)/2),Nt-1)*sdesv;
    intervS=Dmedia+range;
    intervI=Dmedia-range;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box(axes1,'on');
for j=2:size(dados,2)
    cor = [0 0 0];
    corref=[0 0 1];
    if(j==3)
        cor = [0 0 1];
    end
%
    plot(ref(1:jump:end,1),ref(1:jump:end,j),'Parent',axes1,'Color',...
        corref,'LineWidth',3,'MarkerSize',6,'Marker','none',...
        'LineStyle','-','DisplayName','F^{ref}')
%
    plot(ref(1:1:N,1),Dmedia(1:1:end,j-1),'Parent',axes1,'Color',...
        cor,'LineWidth',3,'LineStyle','-','DisplayName','F^{mean}')
%
    if(printa==1)
        plot(ref(1:jumpstat:N,1),intervS(1:jumpstat:end,j-1),'Parent',...
            axes1,'Color',cor,'LineWidth',2,'LineStyle','--')
    %
        plot(ref(1:jumpstat:N,1),intervI(1:jumpstat:end,j-1),'Parent',...
            axes1,'Color',cor,'LineWidth',2,'LineStyle','--')
    end
end
%
Y=yl;
plot(X,Y,'k-')
%
% Create xlabel
xlabel('t','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
xlim(axes1,x1);
% Create ylabel
ylabel('F','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
ylim(axes1,yl);
% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Location','SouthWest','FontSize',12);

% Print
base=['../figuras/prod_' base_name istr];
%print('-djpeg90',base)
print('-depsc','-r300',base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear;