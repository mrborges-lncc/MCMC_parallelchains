clear;
base_name= 'teste';
N=1390;
dados=load('../simuladorBTMM/exp01/conc/conc_ref_0.dat');
file_name = ['../twoStage/select_prod/prod_Max_0.dat'];
dados1=load(file_name);
dados=dados(1:N,:);
dados1=dados1(1:N,:);
%
[num np] = size(dados);
intg=zeros(np-1,1);
intgREF=0.0;
%% calculo da integral da referencia
for i=1:np-1
   intg(i)=trapz(dados(:,1),dados(:,i+1));
   intgREF=intgREF+intg(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------')
disp('INTEGRAL TOTAL DA REFERENCIA')
intgREF
disp('----------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure1 = figure();
%
%% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
   'DataAspectRatio',[max(dados(:,1))*1/1.6 1 2]);
box(axes1,'on');
hold(axes1,'all');

plot(dados(:,1),dados(:,2),'Parent',axes1,'Color',[1 0 0],...
    'MarkerSize',5,'Marker','o','LineStyle','none','DisplayName','ref. 1')
plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],...
   'MarkerSize',5,'Marker','o','LineStyle','none','DisplayName','ref. 2')
plot(dados(:,1),dados(:,4),'Parent',axes1,'Color',[1 1 0],...
   'MarkerSize',5,'Marker','^','LineStyle','none','DisplayName','ref. 3')
plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 0 0],...
   'MarkerSize',5,'Marker','+','LineStyle','none','DisplayName','ref. 4')
plot(dados(:,1),dados(:,6),'Parent',axes1,'Color',[0 1 0],...
   'MarkerSize',5,'Marker','v','LineStyle','none','DisplayName','ref. 5')
plot(dados(:,1),dados(:,7),'Parent',axes1,'Color',[0 1 1],...
   'MarkerSize',5,'Marker','s','LineStyle','none','DisplayName','ref. 6')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num np] = size(dados1);
intga=zeros(np-1,1);
intgAMO=0.0;
% calculo da integral da referencia
for i=1:np-1
   intga(i)=trapz(dados1(:,1),dados1(:,i+1));
   intgAMO=intgAMO+intga(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------')
disp('INTEGRAL TOTAL DA AMOSTRA   ')
intgAMO
disp('----------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculo do error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx=zeros(num,1);
intg2=zeros(np-1,1);
soma=zeros(np-1,1);
for j=1:np-1
   for i=1:num
      aux = (dados(i,j+1)-dados1(i,j+1));
      fx(i)= aux*aux;
   end
%   intg2(j)=trapz(dados1(:,1),fx);
   soma(j)=sum(fx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------')
disp('ERROR                       ')
sum(soma)
%sum(intg2)
disp('----------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');

plot(dados1(:,1),dados1(:,2),'Parent',axes1,'Color',[1 0 0],...
   'LineWidth',2,'DisplayName','data 1')
plot(dados1(:,1),dados1(:,3),'Parent',axes1,'Color',[0 0 1],...
   'LineWidth',2,'DisplayName','data 2')
plot(dados1(:,1),dados1(:,4),'Parent',axes1,'Color',[1 1 0],...
   'LineWidth',2,'DisplayName','data 3')
plot(dados1(:,1),dados1(:,5),'Parent',axes1,'Color',[0 0 0],...
   'LineWidth',2,'DisplayName','data 4')
plot(dados1(:,1),dados1(:,6),'Parent',axes1,'Color',[0 1 0],...
   'LineWidth',2,'DisplayName','data 5')
plot(dados1(:,1),dados1(:,7),'Parent',axes1,'Color',[0 1 1],...
   'LineWidth',2,'DisplayName','data 6')

% Create xlabel
xlabel('t','FontSize',20,'FontName','Times New Roman','FontAngle','italic');%

% Create ylabel
ylabel('C','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',12);

% Print
base=['../figuras/conc_' base_name];
%print('-djpeg90',base)
print('-depsc','-r100',base)

clear
