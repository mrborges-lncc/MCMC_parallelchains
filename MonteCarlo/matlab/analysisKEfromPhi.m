clear all
close all
addpath ~/Dropbox/mrst_borges_tools/
addpath ~/Dropbox/mrst-2021a/
addpath ../../MATLAB_POSPROC/
startup

%% Porosity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = lhsnorm(0,1,50000);
mu_phi  = 0.15;
std_phi = 0.035;  
muY  = log(mu_phi) - (1/2) * log((std_phi^2)/(mu_phi^2) + 1);
varY = log((std_phi^2)/(mu_phi^2) + 1);
stdY = sqrt(varY);
beta = exp(muY);
rho  = sqrt(varY);
p    = beta * exp(rho * Y);
fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('Mean porosity.........: %5.3f\n',mean(p))
fprintf('Std  porosity.........: %5.3f\n',std(p))
fprintf('Mim. porosity.........: %5.3f\n',min(p))
fprintf('Max. porosity.........: %5.3f\n',max(p))
fprintf('=============================================================\n')
fprintf('Rho  value for phi....: %5.3f\n',rho)
fprintf('Beta value for phi....: %5.3f\n',beta)
fprintf('=============================================================\n')
fprintf('=============================================================\n')
%NORMAL(Y,mean(Y),std(Y),'\phi')
%LOGNORMAL(p,mean(p),std(p),'\phi')
%% Permeability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kozeny-Carman
fat = 1/(milli*darcy);
Ss  = 60e03;    % Specific surface area
k_mean = 100/fat;
syms f(x) x
f(x) = x * (mu_phi^3)/(1-mu_phi)^2 - k_mean;
ckc = 1/double(solve(f,x));
c   = 10;            % Kozeny constant
Ss  = sqrt(ckc/c);
k   = (1/(ckc)) * (p.^3) ./ ((1 - p).^2);
mu_k = mean(k);
std_k= std(k);
muY  = log(mu_k) - (1/2) * log((std_k^2)/(mu_k^2) + 1);
varY = log((std_k^2)/(mu_k^2) + 1);
stdY = sqrt(varY);
beta = exp(muY);
rho  = sqrt(varY);
Y    = reverseKlog(k,beta,rho);
k    = beta * exp(rho * Y);
fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('Mean permeability.........: %5.3e (m^2) ... %5.3f (mD)\n',mean(k),mean(k)*fat)
fprintf('Std  permeability.........: %5.3e (m^2) ... %5.3f (mD)\n',std(k),std(k)*fat)
fprintf('Mim. permeability.........: %5.3e (m^2) ... %5.3f (mD)\n',min(k),min(k)*fat)
fprintf('Max. permeability.........: %5.3e (m^2) ... %5.3f (mD)\n',max(k),max(k)*fat)
fprintf('=============================================================\n')
fprintf('Kozeny constant...........: %5.3f\n',c)
fprintf('Specific surface area.....: %5.3f (1/cm)\n',Ss)
fprintf('=============================================================\n')
fprintf('Rho  value for k......: %5.3f\n',rho)
fprintf('Beta value for k......: %5.3e\n',beta)
fprintf('=============================================================\n')
fprintf('=============================================================\n')
plotKC(p,k,1);

%% Young modulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dados obtidos em Flavia (Tese) pg. 134
nu        = 0.20;
c_calcita = 1.302e-11;
c_bulk    = 2.897e-10;
E0_calcita= (1.0/c_calcita) * (3.0*(1.0 - 2*nu));
E         = (1.0/c_bulk) * (3.0*(1.0 - 2*nu));
phi= [0.0:0.05:0.5].';
x  = [0 0.05 0.1 0.15];
mm = [37 25 13];
fm = 2/31;
y  = 0.55*(2e6 + fm * mm * 1e6) * 6894.75729;
E0 = 25.e9;
y  = [E0 y];

syms f a
f(a) = E0 * exp(-a .* x);
lnf(a) = log(f);
syms err(a) derr
err(a)= sum(log(y) - lnf)^2;
derr  =  diff(err,a) == 0;
sp    = double(solve(derr,a));

%p = [0; p];
E = E0 * exp(-sp * p);

Bulk = mean(E)/(3.0*(1.0 - 2*mean(nu)));

mu_E = mean(E);
std_E= std(E);
muY  = log(mu_E) - (1/2) * log((std_E^2)/(mu_E^2) + 1);
varY = log((std_E^2)/(mu_E^2) + 1);
stdY = sqrt(varY);
beta = exp(muY);
rho  = sqrt(varY);

NORMAL(E,mean(E),std(E),'$\mathsf{E}$');

Y = reverseKlog(E,beta,rho);
NORMAL(Y,mean(Y),std(Y),'$Y_{\mathsf{E}}$');
E = beta * exp(rho * Y);

plotEphi(p,E,E0);

fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('Mean Young modulus........: %5.3f (GPa)\n',mean(E)/1e9)
fprintf('Std  Young modulus........: %5.3f (GPa)\n',std(E)/1e9)
fprintf('Mim. Young modulus........: %5.3f (GPa)\n',min(E)/1e9)
fprintf('Max. Young modulus........: %5.3f (GPa)\n',max(E)/1e9)
fprintf('=============================================================\n')
fprintf('E0........................: %5.3f (GPa)\n',E0/1e9)
fprintf('Coefficient a.............: %5.3f \n',sp)
fprintf('=============================================================\n')
fprintf('Rho  value for E......: %5.3f\n',rho)
fprintf('Beta value for E......: %5.3e\n',beta)
fprintf('=============================================================\n')
fprintf('=============================================================\n')



