clear all; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jump = 5;
NC = [1:3];
Ni = [100];
Nf = [200];
Nt = Nf - Ni + 1;
num_datatype = 2;
data_normal  = 0;
expname = 'nCW';
home    = './';
homet   = [home 'data/data'];
homef   = './figuras/';
fileref = ['~/MCMC_parallelchains/twophaseflow/exp/pres/pres_ref_0.dat';
    '~/MCMC_parallelchains/twophaseflow/exp/prod/prod_ref_0.dat'];
%% READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut     = [0 0];
dataref = load_data(fileref,num_datatype,cut);
pres    = dataref{1};
prod    = dataref{2};
X       = [];
vari    = '1';
ch      = 0;
for chain = NC
    ch = ch + 1;
    k  = 0;
    for n = Ni : Nf
        name = [homet expname '_chain' num2str(chain,'%d') '_v' vari...
            '_' num2str(n,'%d') '.dat'];
        aux = load(name,'-mat');
        if k == 0
            sz = size(aux.t);
            Y  = zeros(sz(1),sz(2)-1,Nt);
            X  = aux.t(:,1);
        end
        k = k + 1;
        Y(:,:,k) = aux.t(:,2:end);
    end
end

muY  = mean(Y,3);
stdY = std(Y,0,3);

figdata(pres,muY,stdY,'Pressure ($MPa$)',jump,2)
% Print
base=[homef 'presMCMC_' expname];
%print('-djpeg90',base)
print('-depsc','-r300',base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X    = [];
vari = '2';
ch   = 0;
for chain = NC
    ch = ch + 1;
    k  = 0;
    for n = Ni : Nf
        name = [homet expname '_chain' num2str(chain,'%d') '_v' vari...
            '_' num2str(n,'%d') '.dat']
        aux = load(name,'-mat');
        if k == 0
            sz = size(aux.t);
            Y  = zeros(sz(1),sz(2)-1,Nt);
            X  = aux.t(:,1);
        end
        k = k + 1;
        Y(:,:,k) = aux.t(:,2:end);
    end
end

muY  = mean(Y,3);
stdY = std(Y,0,3);

figdata(prod,muY,stdY,'Oil rate ($m^3/day$)',jump,3)
% Print
base=[homef 'prodMCMC_' expname];
%print('-djpeg90',base)
print('-depsc','-r300',base)

clear all
