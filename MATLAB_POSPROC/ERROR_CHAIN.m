clear;
base_name = 'KL5000_RW_RK0';
%file_name = ['../MCMC/error/erros_' base_name '.dat']
file_name = ['../twoStage/error/erros_' base_name '.dat']
file_name2= ['../twoStage/out/nchain_' base_name '.dat']
dados =load(file_name);
dados2=load(file_name2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[];
y=[];
k=0;
for i=1:size(dados2,1)
    for j=1:dados2(i,2)
        x=[x;k];
        k=k+1;
        y=[y;dados(i,2:end)];
    end
end
dados=[x y];
% Create figure
maximo=0.3; %se max==0 o programa escolhe o valor maximo
maximo2=0.1;
N=size(dados,2);
if(N==3)
    nopt=1;
    A=dados(:,3);
    B=dados(:,3);
end
if(N==4)
    nopt=0;
    A=dados(:,3);
    B=dados(:,4);
end
if(N==5)
    nopt=2;
    A=dados(:,3);
    B=dados(:,4);
end
%figerror2(dados(:,1),dados(:,2),A,B,nopt,maximo)
figerrorColor2(dados(:,1),dados(:,2),A,B,nopt,maximo,maximo2)
base=['../figuras/error_' base_name]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
%print('-djpeg90',base)
size(dados,1)
%clear