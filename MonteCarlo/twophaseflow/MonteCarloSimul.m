clear all
close all

N = 2000;
nome = 'TwoPhase3DMC_only_perm';

currentDir = pwd;
finpe   = [currentDir '/fields/e510x510x20_51x51x5_l50x50x10_'];
finpo   = [currentDir '/fields/e510x510x20_51x51x5_l50x50x10_'];
fperm   = [currentDir '/exp000/fields/perm_amostra_0.dat'];
fporo   = [currentDir '/exp000/fields/poro_amostra_0.dat'];
direct  = dir(fullfile([currentDir '/exp000/fields/*.dat']));

for i = 0:N-1
    n = num2str(i,'%d')
    fper = [finpe num2str(i,'%d') '.dat'];
    fpor = [finpo num2str(i,'%d') '.dat'];
    copyfile(fper, fperm);
    copyfile(fpor, fporo);
    Simulator(0);
    movefile([currentDir '/exp000/conc/sw_amostra_0.dat'],...
        [currentDir '/exp000/conc/sw_' nome '_' n '.dat']);
    movefile([currentDir '/exp000/pres/pres_amostra_0.dat'],...
        [currentDir '/exp000/pres/pres_' nome '_' n '.dat']);
    movefile([currentDir '/exp000/pres/presinj_amostra_0.dat'],...
        [currentDir '/exp000/pres/presinj_' nome '_' n '.dat']);
    movefile([currentDir '/exp000/prod/prod_amostra_0.dat'],...
        [currentDir '/exp000/prod/prod_' nome '_' n '.dat']);
    movefile([currentDir '/exp000/prod/wcut_amostra_0.dat'],...
        [currentDir '/exp000/prod/wcut_' nome '_' n '.dat']);
end