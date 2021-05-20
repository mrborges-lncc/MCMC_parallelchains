clear;
ini=0;
fim=20;
Nthetas=20;
NT =fim-ini+1;
base_name = 'prod_D1_KL20_0';
ref=load('../simul_comp/exp/pressao/pre_ref_0.dat');
%
[FILENAME, PATHNAME] =uigetfile('../twoStage/select_prod/*.dat', 'LOAD DATA');
nc=getNAME(FILENAME)
namet1 = FILENAME(1:nc);
erro=zeros(NT,1);
x=[1:1:NT];
for j=1:NT
    name = [PATHNAME namet1 num2str(ini+j-1) '.dat']
    dados= load(name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    soma=0.0;
    for i=2:size(dados,2)-1
        aux=sum((ref(:,i)-dados(:,i)).^2);
        soma=soma+sqrt(aux);
    end
    erro(j)=soma;
    fprintf('ERRO TOTAL = %f\n',soma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure1 = figure();
B=max(erro);
A=min(erro);
%B=1e-2;
%A=1e-5;
C=min(x);
D=max(x);

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');
plot(x,erro,'Parent',axes1,'Color',[0 0 1],...
    'MarkerSize',4,'LineWidth',2,'DisplayName','erro')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(axes1,[0 D])
ylim(axes1,[A B])
% Create xlabel
xlabel('n','FontSize',20,'FontName','Times New Roman','FontAngle','italic');

% Create ylabel
ylabel('E','FontSize',20,'FontName','Times New Roman','FontAngle','italic');
% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Location','NorthWest','FontSize',12);

% Print
base=['../figuras/conc_' base_name];
%print('-djpeg90',base)
print('-depsc','-r100',base)
clear
