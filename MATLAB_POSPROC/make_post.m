%% DIFFERENTIAL EVOLUTION MARKOV CHAIN algorithm
clear;
close all;
analysis = 10;
tempo = cputime;
d = 4 ;           % dimension of multivariate Gaussian
N = 10000;       % total number of trial
nome = 'DE10';
fat = 16.5;
sig1 = fat*(8.0^2)/21;
sig2 = fat*(8.0^2)/21;
sig3 = fat*(8.0^2)/21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Posterior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cov1 = 0.25*[2 1 0.5 0;
%         1 1.5 0.75 0.5;
%         0.5 0.75 1 0.25;
%         0 0.5 0.25 0.75];
% cov2 = (cov1)^2;
% cov3 = cov1.^(1/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = [-1:0.1:1];
sxd = size(xdata,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a0 = 10;
a1 = -30;
a2 = 10;
a3 = 100;
mu1 = [a3 a2 a1 a0];
a0 = 10;
a1 = 130;
a2 = 10;
a3 = -130;
mu2 = [a3 a2 a1 a0];
a0 = -70.0;
a1 = 20;
a2 = 0;
a3 = -20.0;
mu3 = [a0 a1 a2 a3];
% mu1 = [10 2 30 4];
% mu2 = [-3 3 -5 50];
% mu3 = 10*[5 -4 3 10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1 = [];
theta2 = [];
theta3 = [];
for i=1:N
    Fref1 = [(blackbox(xdata,mu1)+ruido(xdata,sig1))];
    theta1 = [theta1; thetas(xdata,Fref1,d)];
    Fref2 = [(blackbox(xdata,mu2)+ruido(xdata,sig2))];
    theta2 = [theta2; thetas(xdata,Fref2,d)];
    Fref3 = [(blackbox(xdata,mu3)+ruido(xdata,sig3))];
    theta3 = [theta3; thetas(xdata,Fref3,d)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fref1 = [(blackbox(xdata,mu1)+ruido(xdata,0.0))];
Fref2 = [(blackbox(xdata,mu2)+ruido(xdata,0.0))];
Fref3 = [(blackbox(xdata,mu3)+ruido(xdata,0.0))];
tempo = cputime - tempo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref1 = [xdata' Fref1'];
ref2 = [xdata' Fref2'];
ref3 = [xdata' Fref3'];
arqv = '/home/mrborges/MCMC_par/trunk/blackbox/exp/saidaref';
save([arqv '1.dat'],'ref1','-ascii');
save([arqv '2.dat'],'ref2','-ascii');
save([arqv '3.dat'],'ref3','-ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AR=99;
bbfigure12(N,Fref1,Fref2,Fref3,xdata,theta1,theta2,theta3,sxd,AR,mu1,mu2,mu3)
nome = '../figuras/Blackbox_post_';
name = [nome 'curvs'];
print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(analysis==1)
    %N  = 1000;
    cov1 = cov(theta1);
    cov2 = cov(theta2);
    cov3 = cov(theta3);
    m1 = mean(theta1);
    m2 = mean(theta2);
    m3 = mean(theta3);
    v1 = var(theta1);
    v2 = var(theta2);
    v3 = var(theta3);
    theta1 = theta1(end-N+1:end,:,1);
    theta2 = theta2(end-N+1:end,:,1);
    theta3 = theta3(end-N+1:end,:,1);
    nomeh = '../figuras/Blackbox_post_';
    nome = 'curv_'
    name = [nomeh '1'];
    bbfigure(N,mu1,theta1,[nome '1']);
    print('-depsc','-r300',name)
    name = [nomeh '2'];
    bbfigure(N,mu2,theta2,[nome '2']);
    print('-depsc','-r300',name)
    name = [nomeh '3'];
    bbfigure(N,mu3,theta3,[nome '3']);
    print('-depsc','-r300',name)
end
fprintf('\n')