clear;
close all
mu=800.0;
sd=80;
%A=lhsnorm(mu,sd*sd,10000);
Y=lhsnorm(0,1,10000);
A=mu+sd*Y;
MEDIA=mean(A)
VAR=var(A)
hist(A,40);