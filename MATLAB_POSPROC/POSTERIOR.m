close all
clear;
analysis = 1;
tempo = cputime;
d = 4 ;           % dimension of multivariate Gaussian
N = 50000;          % total number of trial
nome = 'DE10';
fat  = 2.0/1.3;
sig1 = (fat*2.5e+00)^2;
sig2 = (fat*2.5e+00)^2;
sig3 = (fat*2.5e+00)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = [-1:0.1:1];
sxd = size(xdata,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Posterior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 = [40 -3 15 12];
mu2 = [20.0 -25 -25 200];
mu3 = [-40 15 -20 60];
% cov1 = 5*covpost(d);
% cov2 = 4*((covpost(d))-2*diag([ 1 0.8 0.6 0.5]'));
% cov3 = 4*sqrt(covpost(d));
cov1 = diag(sig1*ones(d,1));
cov2 = diag(sig2*ones(d,1));%[ 1 0.8 0.6 0.5]');
cov3 = diag(sig3*ones(d,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Criacao dos thetas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1 = mvnrnd(mu1,cov1,N);
theta2 = mvnrnd(mu2,cov2,N);
theta3 = mvnrnd(mu3,cov3,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% verificacao dos thetas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = mean(theta1);
m2 = mean(theta2);
m3 = mean(theta3);
c1 = cov(theta1);
c2 = cov(theta2);
c3 = cov(theta3);
fprintf('\n Estimativa da média de theta1: \n');
m1
fprintf('\n Estimativa da média de theta2: \n');
m2
fprintf('\n Estimativa da média de theta3: \n');
m3
fprintf('\n Estimativa da covariância de theta1: \n');
c1
fprintf('\n Estimativa da covariância de theta2: \n');
c2
fprintf('\n Estimativa da covariância de theta3: \n');
c3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% criacao das curvas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fref1 = [(blackbox(xdata,mu1)+ruido(xdata,0.0))];
Fref2 = [(blackbox(xdata,mu2)+ruido(xdata,0.0))];
Fref3 = [(blackbox(xdata,mu3)+ruido(xdata,0.0))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref1 = [xdata' Fref1'];
ref2 = [xdata' Fref2'];
ref3 = [xdata' Fref3'];
arqv = '/home/mrborges/MCMC_par/trunk/blackbox/exp/saidaref';
save([arqv '1.dat'],'ref1','-ascii');
save([arqv '2.dat'],'ref2','-ascii');
save([arqv '3.dat'],'ref3','-ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AR=100;
home='/home/mrborges/MCMC_par/trunk/'
bbfigure12(N,Fref1,Fref2,Fref3,xdata,theta1,theta2,theta3,sxd,AR,mu1,mu2,mu3)
nome = [home 'figuras/Blackbox_True_post_'];
name = [nome 'curvs'];
print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(analysis==1)
    cov1 = cov(theta1);
    cov2 = cov(theta2);
    cov3 = cov(theta3);
    m1 = mean(theta1);
    m2 = mean(theta2);
    m3 = mean(theta3);
    v1 = var(theta1);
    v2 = var(theta2);
    v3 = var(theta3);
    nome = [ home '/figuras/Blackbox_True_post_'];
    name = [nome '1'];
    bbfigure(N,mu1,theta1,'x');
    print('-depsc','-r300',name)
    name = [nome '2'];
    bbfigure(N,mu2,theta2,'y');
    print('-depsc','-r300',name)
    name = [nome '3'];
    bbfigure(N,mu3,theta3,'z');
    print('-depsc','-r300',name)
end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(analysis==1)
    for i=1:d
        x = theta1(:,i);
        no= ['$x_' num2str(i,3) '$'];
        name = [nome 'x' num2str(i,3)];
        NORMAL(x,mean(x),std(x),no);
        print('-depsc','-r300',name);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:d
        x = theta2(:,i);
        no= ['$y_' num2str(i,3) '$'];
        name = [nome 'y' num2str(i,3)];
        NORMAL(x,mean(x),std(x),no);
        print('-depsc','-r300',name);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:d
        x = theta3(:,i);
        no= ['$z_' num2str(i,3) '$'];
        name = [nome 'z' num2str(i,3)];
        NORMAL(x,mean(x),std(x),no);
        print('-depsc','-r300',name);
    end
end
