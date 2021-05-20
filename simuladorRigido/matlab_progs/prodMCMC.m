clear;
close all;
dados=load('../expref/prod/prod_referencia_0.dat');
%dados=load('../exp/prod/prod_ref_0.dat');
dados=load('../exp/prod/prod_amostra_0.dat');
%dados=load('../prod/prodF_amostra_0.dat');
%dados=load('../../MCMC/select_prod/prodF_6.dat');
%dados=load('../../MCMC/reject_prod/prodF_0.dat');
nprod=size(dados,2)-1;
%
for i=2:nprod+1
    name = ['well ' num2str(i-1)]
    fig_conc(dados(:,1),dados(:,i),name)
end
clear