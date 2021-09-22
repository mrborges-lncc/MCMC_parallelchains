clear all
close all

Nch_ini = 0;
Nch_fim = 5;
Nchains = (Nch_fim - Nch_ini) + 1;
Nini = repmat(0, 1, Nchains);
Nini = [1255 1422 1377 1366 1334 1272];
Nfim = repmat(0, 1, Nchains);
Nfim = [1363 1547 1522 1498 1472 1417];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];

nome= 'TwoPhase3D_RW_RK';
%% Parallel
Npar = 2;
cont = 0;
while mod(Nchains,Npar) ~= 0 && cont < 100
    cont = cont + 1;
    Npar = Npar - 1;
end

fprintf('\n Number of process: %d\n\n',Npar);

delete(gcp('nocreate'));
parpool('local',Npar);

currentDir = pwd;

%for k = Nch_ini : Nch_fim
parfor (j = 0 : Nchains-1, Npar)
    k = j;
    m = k + 1;
    rk= num2str(10+k,'%4.3d');
    fprintf('\nExp no. %d %s\n',k,rk);
    finpe = [currentDir(1:end-23) 'twoStage/select_fields/field_v1_' nome num2str(k,'%d') '_'];
    finpo = [currentDir(1:end-23) 'twoStage/select_fields/field_v2_' nome num2str(k,'%d') '_'];
    for i = Nini(m):1:Nfim(m)
        n = num2str(i,'%d');
        fper = [finpe num2str(i,'%d') '.dat'];
        fpor = [finpo num2str(i,'%d') '.dat'];
        fperm= [currentDir '/exp' rk '/fields/perm_amostra_0.dat'];
        fporo= [currentDir '/exp' rk '/fields/poro_amostra_0.dat'];
        copyfile(fper, fperm);
        copyfile(fpor, fporo);
        SimulatorMCMC(str2num(rk));
        movefile([currentDir '/exp' rk '/conc/sw_amostra_0.dat'],...
            [currentDir '/exp/conc/sw_MCMC_' nome num2str(k,'%d') '_' n '.dat']);
        movefile([currentDir '/exp' rk '/pres/pres_amostra_0.dat'],...
            [currentDir '/exp/pres/pres_MCMC_' nome num2str(k,'%d') '_' n '.dat']);
        movefile([currentDir '/exp' rk '/pres/presinj_amostra_0.dat'],...
            [currentDir '/exp/pres/presinj_MCMC_' nome num2str(k,'%d') '_' n '.dat']);
        movefile([currentDir '/exp' rk '/prod/prod_amostra_0.dat'],...
            [currentDir '/exp/prod/prod_MCMC_' nome num2str(k,'%d') '_' n '.dat']);
        movefile([currentDir '/exp' rk '/prod/wcut_amostra_0.dat'],...
            [currentDir '/exp/prod/wcut_MCMC_' nome num2str(k,'%d') '_' n '.dat']);
    end
end
delete(gcp('nocreate'));
clear