clear;
base_name = 'KL20';
%file_name = ['../MCMC/error/erros_' base_name '.dat']
file_name = ['../twoStage/out/lambda_' base_name '.dat']
dados=load(file_name);
% Create figure
maximo=3000; %se max==0 o programa escolhe o valor maximo
minimo=1000;
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
figlambdas(dados(:,1),dados(:,2),A,B,nopt,maximo,minimo);
base=['../figuras/lambda_' base_name];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
%print('-djpeg90',base)
size(dados,1)
clear