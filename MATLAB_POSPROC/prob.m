close all
clear;
analysis = 1;
tempo = cputime;
d = 11;           % dimension of multivariate Gaussian
N = 5000;          % total number of trial
nome = 'DE10';
fat  = 2.0e-02/1.3;
sig1 = 1;
sig2 = 1;
sig3 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = [-2:0.1:2];
sxd = size(xdata,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Posterior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 = [0.14 7e-15 -1.3 -5.3e-14 3.6 1.3e-13 -2.7 -1.2e-13 -1.3 2.9e-14 0.53];
mu2 = [-9.6e-02 0 0.75 0 -1.4 0 -1 0 2.5 0 3];
mu3 = -mu1;
mu3(end)=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% criacao das curvas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fref1 = [(blackbox(xdata,mu1)+ruido(xdata,0.0))];
Fref2 = [(blackbox(xdata,mu2)+ruido(xdata,0.0))];
Fref3 = [(blackbox(xdata,mu3)+ruido(xdata,0.0))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = zeros(N,d);
x2 = zeros(N,d);
x3 = zeros(N,d);
for i=1:N
    y1 = Fref1+sig1*randn(1,size(Fref1,2));
    x1(i,:) = regres10(xdata,y1);
    y2 = Fref2+sig2*randn(1,size(Fref2,2));
    x2(i,:) = regres10(xdata,y2);
    y3 = Fref3+sig3*randn(1,size(Fref3,2));
    x3(i,:) = regres10(xdata,y3);
end
sample1 = zeros(N,size(Fref1,2));
sample2 = zeros(N,size(Fref2,2));
sample3 = zeros(N,size(Fref3,2));
for i=1:N
    sample1(i,:) = blackbox(xdata,x1(i,:));
    sample2(i,:) = blackbox(xdata,x2(i,:));
    sample3(i,:) = blackbox(xdata,x3(i,:));
end
sigma1 = std(sample1);
sigma2 = std(sample2);
sigma3 = std(sample3);
plot(xdata,Fref1,'-k','LineWidth',3)
hold on
plot(xdata,Fref2,'-r','LineWidth',3)
plot(xdata,Fref3,'-b','LineWidth',3)
errorbar(xdata,blackbox(xdata,mean(x1)),sigma1,'Marker','o','MarkerSize',8,...
    'LineWidth',2,'Color',[0 0 0],'LineStyle','none');
errorbar(xdata,blackbox(xdata,mean(x2)),sigma2,'Marker','s','MarkerSize',8,...
    'LineWidth',2,'Color',[1 0 0],'LineStyle','none');
errorbar(xdata,blackbox(xdata,mean(x3)),sigma3,'Marker','^','MarkerSize',8,...
    'LineWidth',2,'Color',[0 0 1],'LineStyle','none');
