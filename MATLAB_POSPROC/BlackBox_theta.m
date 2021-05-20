%% DIFFERENTIAL EVOLUTION MARKOV CHAIN algorithm
clear;
close all;
analysis = 1;
tempo = cputime;
d = 4 ;% dimension of multivariate Gausesian
Nc0=0;
Nc= 7;            % Number of chains
NT=Nc-Nc0+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nini=00;
Nfim=1;
N=(Nfim-Nini)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
tm   = 0;
base_name = 'blackbox_RK';
for i=Nc0:Nc
    file_name = ...
        ['../twoStage/error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
end
tm=tm-1;
%N=max(tm);
clear dados
ref1=load('../blackbox/exp/saidaref1.dat');
ref2=load('../blackbox/exp/saidaref2.dat');
ref3=load('../blackbox/exp/saidaref3.dat');
variavel = '../twoStage/select_thetas/theta_v';
variavel1 = '../twoStage/out/nchain_';
X1=zeros(N,d,Nc);
X2=zeros(N,d,Nc);
X3=zeros(N,d,Nc);
%X=zeros(M,d,Nc);
nchain=zeros(N,Nc);
nk=0;
SOMA = [];
for i=Nc0:Nc
    nk = 0;
    nfile =...
        [variavel1 base_name num2str(i,3) '.dat'];
    d = load(nfile);
    SOMA = [SOMA; sum(d(:,2))];
    nchain(1:size(d,1),i+1) = d(:,2);
    for k=Nini:Nfim
        nfile1 =...
            [variavel '1_' base_name num2str(i,3) '_' num2str(k,5) '.dat'];
        nfile2 =...
            [variavel '2_' base_name num2str(i,3) '_' num2str(k,5) '.dat'];
        nfile3 =...
            [variavel '3_' base_name num2str(i,3) '_' num2str(k,5) '.dat'];
        nk=nk+1;
        X1(nk,:,i+1) = load(nfile1)';
        X2(nk,:,i+1) = load(nfile2)';
        X3(nk,:,i+1) = load(nfile3)';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loc = 1;
% X=[];
% for i=0:Nc-1
%     loc = 1;
%     for j=1:tm(i+2)
%         loc1 = loc;
%         loc2 = loc+nchain(j,i+1)-1;
%         for k = loc1:loc2
%             X = [X ; X1(j,:,i+1)];
%         end
%         loc = loc2+1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = [-1:0.1:1];
sxd = size(xdata,2);
a0 = 10;
a1 = 3;
a2 = -5;
a3 = -1.25;
mu1 = [a0 a1 a2 a3];
mu2 = [a3 a2 a1 a0];
a0 = 150;
a1 = 25;
a2 = -25;
a3 = 80;
mu3 = [a0 a1 a2 a3];
sig = 0.0001;
Fref1 = [(blackbox(xdata,mu1)+ruido(xdata,sig))];
Fref2 = [(blackbox(xdata,mu2)+ruido(xdata,sig))];
Fref3 = [(blackbox(xdata,mu3)+ruido(xdata,sig))];
sigma = 0.05;
sigma1 = 0.045;
sigma2 = 0.045;
sigma3 = 0.045;
%ss1 = eye(d,d); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempo = cputime - tempo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bbfigure2(Er)
% nome = './figuras/blackbox3_';
% name = [nome 'error'];
% print('-depsc','-r300',name)

bbfigure12(N,Fref1,Fref2,Fref3,xdata,X1,X2,X3,sxd)
nome = '../figuras/blackbox_';
name = [nome 'curvs'];
print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(analysis==1)
    %N  = 1000;
    X1 = X1(end-N+1:end,:,Nc+1);
    X2 = X2(end-N+1:end,:,Nc+1);
    X3 = X3(end-N+1:end,:,Nc+1);
    nome = '../figuras/blackbox_';
    name = [nome '1'];
    bbfigure(N,mu1,X1);
    print('-depsc','-r300',name)
    name = [nome '2'];
    bbfigure(N,mu2,X2);
    print('-depsc','-r300',name)
    name = [nome '3'];
    bbfigure(N,mu3,X3);
    print('-depsc','-r300',name)
end
% sz=size(Er,2);
% sz=sz-ceil(sz*0.25);
% Y = Er(sz:end);
% mu= mean(Er(sz:end));
% std=std(Er(sz:end));
% tipo = 'Er'
% NORMAL(Y,mu,std,tipo);
% nome = '../figuras/blackbox_';
% name = [nome 'hist'];
% print('-depsc','-r300',name)
fprintf('\n################################################################\n')
%Taxa = (cont/double(N))*100.0
fprintf('\n')