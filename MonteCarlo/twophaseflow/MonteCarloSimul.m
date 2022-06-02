clear all
close all

delete(gcp('nocreate'))

ini = 0;
fim = 1;
N   = fim - ini + 1;
nome = 'TwoPhase3D_KC_MC';

home = '~/fields/campos/';

currentDir = home;
finpe   = [currentDir 'perm_'];
finpo   = [currentDir 'phi_'];
finyo   = [currentDir 'E_'];
currentDir = pwd;
direct  = dir(fullfile([currentDir '/exp/fields/*.dat']));

Npar = 4;

parpool('local',Npar);
for i = ini:Npar:fim
    parfor (j = 0:Npar-1, Npar)
%    for j = 0:Npar-1
        n = j;
        k = i + j
        fper = [finpe num2str(k,'%d') '.dat'];
        fpor = [finpo num2str(k,'%d') '.dat'];
        fyou = [finyo num2str(k,'%d') '.dat'];
        fperm= [currentDir '/exp' num2str(n,'%4.3d') '/fields/perm_amostra_0.dat'];
        fporo= [currentDir '/exp' num2str(n,'%4.3d') '/fields/poro_amostra_0.dat'];
        fyoun= [currentDir '/exp' num2str(n,'%4.3d') '/fields/young_amostra_0.dat'];
        
        fprintf('\n%d %d\n',n,k)
        fprintf('%s\n',fyou)
        fprintf('%s\n',fperm)
        fprintf('%s\n',fporo)
        fprintf('%s\n',fyoun)
        
        copyfile(fper, fperm);
        copyfile(fpor, fporo);
        copyfile(fyou, fyoun);
        fprintf('%s => %s\n',fyou,fyoun)
        Simulator(n);
        movefile([currentDir '/exp' num2str(n,'%4.3d') '/conc/sw_amostra_0.dat'],...
            [currentDir '/exp/conc/sw_' nome '_' num2str(k,'%d') '.dat']);
        movefile([currentDir '/exp' num2str(n,'%4.3d') '/pres/pres_amostra_0.dat'],...
            [currentDir '/exp/pres/pres_' nome '_' num2str(k,'%d') '.dat']);
        movefile([currentDir '/exp' num2str(n,'%4.3d') '/pres/presinj_amostra_0.dat'],...
            [currentDir '/exp/pres/presinj_' nome '_' num2str(k,'%d') '.dat']);
        movefile([currentDir '/exp' num2str(n,'%4.3d') '/prod/prod_amostra_0.dat'],...
            [currentDir '/exp/prod/prod_' nome '_' num2str(k,'%d') '.dat']);
        movefile([currentDir '/exp' num2str(n,'%4.3d') '/prod/wcut_amostra_0.dat'],...
            [currentDir '/exp/prod/wcut_' nome '_' num2str(k,'%d') '.dat']);
        fprintf('\n\n %s\n',[currentDir '/exp' num2str(n,'%4.3d') '/conc/sw_amostra_0.dat'])
        fprintf('\n\n %s\n',[currentDir '/exp/conc/sw_' nome '_' num2str(k,'%d') '.dat'])
    end
end
delete(gcp('nocreate'))