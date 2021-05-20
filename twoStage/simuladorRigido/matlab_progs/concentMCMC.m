clear;
close all
dados=load('../expref/conc/conc_referencia_0.dat');
%dados=load('../exp/conc/conc_ref_0.dat');
% dados=load('../exp/conc/conc_amostra_0.dat');
nprod=size(dados,2)-1;
%
for i=2:nprod+1
    name = ['p' num2str(i-1)]
    fig_conc(dados(:,1),dados(:,i),name)
end
clear