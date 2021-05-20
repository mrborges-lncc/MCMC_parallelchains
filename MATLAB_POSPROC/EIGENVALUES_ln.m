clear;
base_name = 'KL';
M=10000;
base_name = [base_name num2str(M,5)]
file_name1 = ['../gera_KL/MATLAB/out/aval1_1x1_100x100_l0.2_M10000.dat'];
file_name2 = ['../gera_KL/MATLAB/out/aval3_1x1_100x100_l0.2_M10000.dat'];
file_name = ['../gera_KL/MATLAB/out/aval2_1x1_100x100_b0.5_M10000.dat'];
dados1=load(file_name1);
dados2=load(file_name2);
total=sum(dados1(1:end));
trun =sum(dados1(1:M));
strength = (trun/total);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = ['Energia obtida com ' num2str(M,'%8.2d') ' termos na expansao = '...
    num2str(100*strength,'%8.2f')];
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(nome)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figavalln(M,dados1,dados2,strength);
base=['../figuras/lnaval_' base_name];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
%print('-djpeg90',base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo do strength
lim_st = 0.98;
for i=1:size(dados1)
    total=sum(dados1(1:end));
    trun =sum(dados1(1:i));
    strength = (trun/total);
    if(strength>=lim_st)
        M=i;
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = ['Numero de termos na expansao KL (M) para obter uma energia de '...
    num2str(strength,'%8.2f') ' = ' num2str(M,'%8.2d')];
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(nome)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
clear