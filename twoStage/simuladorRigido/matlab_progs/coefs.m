rho=1.0;
my=0.0;
vy=1.0;
M=8e-14;
Y=randn(1,1000);
disp('############################')
disp('####----- CAMPO Y ------####')
%my=mean(Y)
%vy=var(Y)
maxy=max(Y);
miny=min(Y);
%hist(Y,50)
%%%%%%%%%%%%%%%%%%%%
mE=exp(rho*my+rho*rho*vy/2);
m=M/mE;
mE=m*exp(rho*my+rho*rho*vy/2);
vE=m*m*exp(2*rho*my+rho*rho*vy)*(exp(rho*rho*vy)-1);
stdE=sqrt(vE);
E=m*exp(rho*Y);
disp('############################')
disp('####----- CAMPO E ------####')
maxE=max(E)
minE=min(E)
mE=mean(E)
vE=var(E)
stdE=sqrt(vE)
%hist(E,50)
disp('############################')
disp('####----- VALORES ------####')
m
rho
disp('############################')
clear
