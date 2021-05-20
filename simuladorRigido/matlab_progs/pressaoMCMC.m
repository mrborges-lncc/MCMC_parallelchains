clear;
dados=load('../expref/prod/prod_referencia_0.dat');
dados=load('../exp/pres/pres_ref_0.dat');
dados=load('../exp/pres/pres_amostra_0.dat');
nprod=size(dados,2)-1;
%
for i=2:nprod+1
    name = ['p' num2str(i-1)]
    fig_pres(dados(2:end,1),dados(2:end,i),name)
end
clear