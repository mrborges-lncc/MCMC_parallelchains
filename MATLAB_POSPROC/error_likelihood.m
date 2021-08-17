clear;
close all
ini = 11;
fim = 75;
sig = 1e-03;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = 'TwoPhase3D_onlyPerm_RW_RK';
nome = 'TwoPhase3D_RW_RK';
base_name = ['prod_D2_' nome];
hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
% hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMCrw_onlyPerm/';
hom = '../';
homf= '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'twophaseflow/exp/prod/prod_ref_0.dat']);
% dados=load([hom 'twophaseflow/exp/pres/pres_ref_0.dat']);
ref=dados(ini:fim,2:end);
% dados=load([hom 'twophaseflow/exp000/prod/prod_amostra_0.dat']);
dados=load([hom 'twoStage/select_prod/' base_name '0_0.dat']);
dat=dados(ini:fim,2:end);
norma_ref = (norm(ref).^2)
norma_l2  = sum(sum((ref - dat).^2))
erro = (norma_l2)/norma_ref
lratio = (0.5/sig)*sqrt(erro)
lratio = (1/(sig))*erro

clear dados
