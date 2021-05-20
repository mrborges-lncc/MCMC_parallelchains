clear;
close all
%dados=load('../exp/conc/conc_referencia_0.dat');
dados=load('../exp/conc/conc_ref_0.dat');
%dados1=load('../exp/conc/conc_amostra_0.dat');
%dados1=load('../exp/conc/conc_ref_0.dat');
dados1=load('../../twoStage/select_prod/prod_D1_FS_DE_RK3_709.dat');
%dados=load('../../forecast/simuladorRigido/exp/prod/prod_MCMC_2.dat');
%dados1=load('../../forecast/simuladorRigido/exp/conc/conc_MC_2.dat');
nprod=size(dados,2)-1;
NN = size(dados,1)-1;
%
soma = 0.0;
denom = 0.0;
for i=2:size(dados,2)
    soma = soma + norm(dados(1:NN,i)-dados1(1:NN,i),'inf')/norm(dados(1:NN,i),'inf');
end
NORMA=soma
%
for i=2:nprod+1
    name = ['p' num2str(i-1)];
    fig_concompara(dados(1:NN,1),dados(1:NN,i),...
        dados1(1:NN,1),dados1(1:NN,i),name)
end
%clear