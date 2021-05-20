clear;
close all
%dados=load('../expref/conc/conc_referencia_0.dat');
dados=load('../exp/prod/prod_ref_0.dat');
%dados1=load('../exp/prod/prod_ref_0.dat');
%dados1=load('../exp/prod/prod_amostra_0.dat');
dados1=load('../../twoStage/select_prod/prod_D1_FS_DREAM_RK1_2.dat');
nprod=size(dados,2)-1;
%
for i=2:nprod+1
    name = ['p' num2str(i-1)]
    fig_concompara(dados(:,1),dados(:,i),...
        dados1(:,1),dados1(:,i),name)
end
clear