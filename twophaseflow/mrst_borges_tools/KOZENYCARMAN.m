%--------------------------------------------------------------------------
%
% PROGRAMA PARA A CRIACAO DOS CAMPOS GAUSSIANOS DA PERMEABILIDADE A PARTIR%
% DO CAMPO DE POROSIDADES 
%
close all;%
clear all;%
TOL     = 0.05;
TOLZERO = 0.0075;
corretor=1;
N=1;
stat=1;
nome= 'stat10amostras';
nome= 'statREF';
home = '~/projMCMC/trunk/figuras/';
%home = '~/GEO_MCMC/trunk/figuras/';
home = '/home/mrborges/MCMC_par/trunk/SIMULADOR_VISCOELASTICO/figuras/'
%home = '/home/mrborges/MCMC_exp/trunk/simuladorRigido/figuras/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIPO DE RELACAO ENTRE PHI e E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
self = 10; % se igual usa relacao self consistent | c.c. relacao de Sprigg
tipophi = 10; % if 1 \phi gaussian
             % else Coefficients of Porosity field => \phi = M*\exp[\alpha*Y]
phinormalpadrao = 1; % if 1 normal padrao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEDIAPHI=0.20;
STDPHI  =sqrt(0.0010);%0.0866;
MAXPHI  =.75;
MINPHI  =0.005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(self==1)
    beta = 0.9875;
    stdKs= 0.0e00;
    medE = 5.55e09
    stdE = 0.00e08;
    nu   = 0.13;
    Kmed = medE/(3*(1-2*nu));
%    Kmed = 5.0e09;
    Ks   = Kmed/(1.0-beta)
    G    =  (3*Ks*(1-2*nu))/(2*(1+nu));
    Kfw  = 0.0*2.2e09;
    Kfo  = 0.0*2.2e09;
    Kf   = (Kfw+Kfo)/2;
else
%     SE    = 3.475e09;
%     stdSE = 0.05e09;
%     sY    = 8.5;
    SE    = 1.90e10;
    stdSE = 0.1e10;
    stdSE = 0.0e10;
    sY    = 4.0;
    SE    = 1.1e10;
    stdSE = 0.0025e10;
    ME    = 1.152e09;
    SE    = 1.1520e11;
%    SE    = 7.68e10
    stdSE = 0.2e10;
    %ME    = 1.0e10;
    sY    = (log(SE)-log(ME))/MEDIAPHI
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific surface [m^-1]
KMED = 9.869233e-14;
%KMED = 8.55e-14;
%KMED = 1.0
ss=6.8e5;
ss=4.48e5;
ss=3.75e05;
%ss=sqrt((MEDIAPHI^3)/((1-MEDIAPHI)^2))*0.99995
co=1;
M=50000;
for i=1:M
    ss=ss+ss*(i*1e-7);
    k= sum(kcm(co,ss,MEDIAPHI+STDPHI*randn(M,1)))/M
    k=(k - KMED)/KMED;
    if(abs(k)<5e-4)
        SS=ss;
        break
    end
end
%SS = 4.828e05;
%SS = 2.675e05;
%SS = 2.3695e+05
SS = 3.75e+05
stdSS = SS*0.025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VARPHI=STDPHI^2;
muY   =log(MEDIAPHI)-0.5*(log((VARPHI/(MEDIAPHI*MEDIAPHI))+1));
varY  =log((VARPHI/(MEDIAPHI*MEDIAPHI))+1);
alpha =sqrt(varY);
M     =exp(muY);
muX   =exp(muY+varY/2);
varX  =exp(2*muY+varY)*(exp(varY)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux  = '../figuras/';
base_aux2 = '/home/mrborges/MCMC_par/trunk/SIMULADOR_VISCOELASTICO/expref/fields/';
base_aux2 = '/home/mrborges/MCMC_exp/trunk/simuladorRigido/exp/fields/';
% base_aux2 = '/home/mrborges/projMCMC/trunk/gera_KL/MATLAB/campos/';
%base_aux2 = '/home/mrborges/projMCMC/trunk/SIMULADOR_VISCOELASTICO/exp/fields/';
%base_aux2 = '/home/mrborges/GEOMECH/TwoWay_0/trunk/SIMULADOR_VISCOELASTICO/exp/fields/';
% base_aux2 = '/home/mrborges/teste/GEO_MCMC/trunk/SIMULADOR_VISCOELASTICO/exp/fields/';
% base_aux2 = '/mnt/dados/mrborges/GEO_MCMC/trunk/SIMULADOR_VISCOELASTICO/exp/fields/';
base_aux2 = '~/TwoWay_lastversionPAPER/trunk/SIMULADOR_VISCOELASTICO/expSISMELAST/fields/';
imprK=10; % se igual a 1 gera o grafico de K

disp('----------------')
disp(' LOADING FIELD ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/Dropbox/KLE/fields/e*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('/home/mrborges/TwoWay_lastversionPAPER/trunk/SIMULADOR_VISCOELASTICO/expSISMELAST/fields/phi_ref*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile(...
%'/home/mrborges/projMCMC/trunk/gera_KL/MATLAB/camposphi/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile(...
%'/home/mrborges/GEOMECH/TwoWay_0/trunk/SIMULADOR_VISCOELASTICO/exp/fields/*.dat', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file=/home/mrborges/Mcoll/field_klfrac/Yfrac05_s1_M100_120x120_*.dat
n=0;
for i=size(line_file,2):-1:1
    a=line_file(i);
    if(a=='.')
        in=line_file(i-1);
        for j=i-2:-1:1
            a=line_file(j);
            if(a=='_')
                str_k=j;
                i=0;
                break
            else
                in=[a in];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LEITURA DOS CAMPOS GAUSSIANOS DE PHI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_base=line_file(1:str_k);
MAX=-1e30;
MIN=1e30;
MEDIAK=0.0;
MMEDIO=0.0;
ALPHAMEDIA=0.0;
MEDIAE=0.0;
MEDIAP=0.0;
varTotalG=[];
varTotalP=[];
varTotalK=[];
varTotalE=[];
ini=str2num(in);
logmodY=[];
logK=[];
permeabilidade = [];
porosidade = [];
young = [];
vetMM = [];
veteM = [];
vetaalpha= [];
vetealpha= [];
fim=ini+N-1;
for II=ini:1:fim
    line_file = [file_base num2str(II,5) '.dat']
    fid = fopen(line_file,'r');
    mattamp = fscanf(fid,'%f');

    disp('file loaded.')
    fclose(fid);
    inf = mattamp(1:8);
    Lx = inf(1);
    Ly = inf(2);
    nx = inf(3);
    ny = inf(4);
    ntipo = inf(5);
    beta = inf(6);
    dx = Lx/nx;
    dy = Ly/ny;
    if(abs(nx-ny)<1e-8)
        quad=1;
    else
        quad=0;
    end
    mattamp = mattamp(9:length(mattamp));
    permmap=[];
    k=0;
    for j=ny:-1:1
        k=k+1;
        if(mattamp(k)~=ny-j)
            disp('erro1')
            break
        end
        for i=1:nx
            k=k+1;
            permmap(j,i)=mattamp(k);
        end
        k=k+1;
        if(mattamp(k)~=192837465)
            disp('erro2')
            break
        end       
    end
    clear mattamp
    s = size(permmap);
    x = s(:,2);
    y = s(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CORRECAO DA MEDIA E VARIANCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd);
    media=mean(vd);
    std=sqrt(variancia);
    if(corretor==1)
        vd=reshape(permmap,nx*ny,1);
        variancia=var(vd);
        permmap=sqrt(1.0/variancia)*permmap;
        vd=reshape(permmap,nx*ny,1);
        media=mean(vd);
        permmap=(-media)+permmap;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SALVANDO O CAMPO COM MEDIA ZERO E VARIANCIA UM %%%%%%%%%%%%%%%
    base1=['phi_' base_name];
    nfile= num2str(II,5);
    bb2 = [base_aux2 base1 '_' nfile '.dat']
%         bb2 = line_file;
    outfile2 = fopen(bb2, 'wt');
    fprintf(outfile2,'%f\n',Lx);
    fprintf(outfile2,'%f\n',Ly);
    fprintf(outfile2,'%d\n',nx);
    fprintf(outfile2,'%d\n',ny);
    fprintf(outfile2,'%d\n',ntipo);
    fprintf(outfile2,'%f\n',beta);
    fprintf(outfile2,'%d\n',2);
    fprintf(outfile2,'%d\n',2);
%
    cont=0
    for j=ny:-1:1
        fprintf(outfile2,'%d\n',ny-j);
        for i=1:nx
            if(tipophi == 1)
                yy = STDPHI*permmap(j,i)+MEDIAPHI;
                if(yy<TOLZERO)
                    disp('*************************************');
                    disp('Porosidade menor que zero (ou limite)');
                    disp('*************************************');
                    yy = TOLZERO;
                    cont=cont+1;
                    permmap(j,i) = (yy-MEDIAPHI)/STDPHI;
                end
                permmap(j,i) = yy;
            end
            fprintf(outfile2,'%f ',permmap(j,i));
        end
      fprintf(outfile2,'\n192837465\n');
    end
    fclose(outfile2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd);
    media=mean(vd);
    std=sqrt(variancia);
    varTotalG=[varTotalG; vd];
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DO CAMPO %d GAUSSIANO DA POROSIDADE\n',II)
    fprintf('MEDIA...........: %f\n',media);
    fprintf('VARIANCIA.......: %f\n',variancia);
    fprintf('DESVIO PADRAO...: %f\n',std);
    fprintf('MAXIMO..........: %f\n',max(max(permmap)));
    fprintf('MINIMO..........: %f\n',min(min(permmap)));
    fprintf('#########################################################\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(tipophi~=1)
        permmap=M*exp(alpha*permmap);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vdp=reshape(permmap,nx*ny,1);
    media=mean(vdp);
    MEDIAP=MEDIAP+media;
    MAX=max(MAX,max(vdp));
    MIN=min(MIN,min(vdp));
    if(stat==1)
        porosidade = [porosidade; vdp];
    end
    variancia=var(vdp);
    std=sqrt(variancia);
    varTotalP=[varTotalP; vdp];
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DO CAMPO %d DE POROSIDADES\n',II)
    fprintf('MEDIA...........: %f\n',media);
    fprintf('VARIANCIA.......: %f\n',variancia);
    fprintf('DESVIO PADRAO...: %f\n',std);
    fprintf('MAXIMO..........: %f\n',max(max(permmap)));
    fprintf('MINIMO..........: %f\n',min(min(permmap)));
    fprintf('#########################################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CRIACAO DE K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perm=[];
    stdK = 1.e-15;
    for j=1:ny
        for i=1:nx
            perm(j,i) = kcm(co,SS+stdSS*randn(1,1),permmap(j,i));
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vdk=reshape(perm,nx*ny,1);
    media=mean(mean(vdk));
    MEDIAK=MEDIAK+media;
    if(stat==1)
        permeabilidade=[permeabilidade; vdk];
    end
    variancia=var(vdk);
    std=sqrt(variancia);
    varTotalK=[varTotalK; vdk];
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DO CAMPO %d DE PERMEABILIDADES\n',II)
    fprintf('MEDIA...........: %e\n',media);
    fprintf('VARIANCIA.......: %e\n',variancia);
    fprintf('DESVIO PADRAO...: %e\n',std);
    fprintf('MAXIMO..........: %e\n',max(max(perm)));
    fprintf('MINIMO..........: %e\n',min(min(perm)));
    fprintf('#########################################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    base2=['perm_' base_name];
%******* impressao dos campos LN(K) ***************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******* impressao dos campos K(PHI) ** COZENY CARMAN *********************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    muY   =log(media)-0.5*(log((variancia/(media*media))+1));
    varY  =log((variancia/(media*media))+1);
    aalpha =sqrt(varY);
    vetaalpha = [vetaalpha; aalpha];
    MM     = exp(muY);
    vetMM  = [vetMM; MM];
    nfile= num2str(II,5);
    bb2 = [base_aux2 base2 '_' nfile '.dat']
    outfile2 = fopen(bb2, 'wt');
    fprintf(outfile2,'%f\n',Lx);
    fprintf(outfile2,'%f\n',Ly);
    fprintf(outfile2,'%d\n',nx);
    fprintf(outfile2,'%d\n',ny);
    fprintf(outfile2,'%d\n',ntipo);
    fprintf(outfile2,'%f\n',beta);
    fprintf(outfile2,'%d\n',2);
    fprintf(outfile2,'%d\n',2);
%
    for j=ny:-1:1
        fprintf(outfile2,'%d\n',ny-j);
        for i=1:nx
            yy = (log(perm(j,i))-log(MM))/aalpha;
            if(stat==1)
                logK = [logK ; yy];
            end
            fprintf(outfile2,'%f ',yy);
        end
        fprintf(outfile2,'\n192837465\n');
    end
    fclose(outfile2);
    clear perm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CRIACAO DE E (modulo de Young) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E=[];
    for j=1:ny
        for i=1:nx
            if(self==1)
                E(j,i) = selfcons(permmap(j,i),G,Ks,Kf,nu) + stdE*randn(1,1);
            else
                 E(j,i)=sc(SE + stdSE*randn(1,1),sY,permmap(j,i));
            end
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vdE=reshape(E,nx*ny,1);
    media=mean(vdE);
    MEDIAE=MEDIAE+media;
    if(stat==1)
        young = [young; vdE];
    end
    variancia=var(vdE);
    std=sqrt(variancia);
    varTotalE=[varTotalE; vdE];
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DO CAMPO %d MODULO DE YOUNG\n',II)
    fprintf('MEDIA...........: %e\n',media);
    fprintf('VARIANCIA.......: %e\n',variancia);
    fprintf('DESVIO PADRAO...: %e\n',std);
    fprintf('MAXIMO..........: %e\n',max(max(E)));
    fprintf('MINIMO..........: %e\n',min(min(E)));
    fprintf('#########################################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    base2=['young_' base_name];
%******* impressao dos campos LN(E) ***************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******* impressao dos campos E(PHI) ** SELF CONSISTENCY ******************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        muY   =log(media)-0.5*(log((variancia/(media*media))+1));
        varY  =log((variancia/(media*media))+1);
        ealpha =sqrt(varY);
        vetealpha = [vetealpha; ealpha];
        eM     =exp(muY);
        veteM  = [veteM; eM];
        nfile= num2str(II,5);
        bb2 = [base_aux2 base2 '_' nfile '.dat']
        outfile2 = fopen(bb2, 'wt');
        fprintf(outfile2,'%f\n',Lx);
        fprintf(outfile2,'%f\n',Ly);
        fprintf(outfile2,'%d\n',nx);
        fprintf(outfile2,'%d\n',ny);
        fprintf(outfile2,'%d\n',ntipo);
        fprintf(outfile2,'%f\n',beta);
        fprintf(outfile2,'%d\n',2);
        fprintf(outfile2,'%d\n',2);
%
        for j=ny:-1:1
          fprintf(outfile2,'%d\n',ny-j);
          for i=1:nx
              yy = (log(E(j,i))-log(eM))/ealpha;
              fprintf(outfile2,'%f ',yy);
              if(stat==1)
                  logmodY = [logmodY; yy];
              end
          end
          fprintf(outfile2,'\n192837465\n');
        end
        fclose(outfile2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(stat==1)
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DOs CAMPOs log-MODULO DE YOUNG\n')
    fprintf('MEDIA...........: %e\n',mean(logmodY));
    fprintf('VARIANCIA.......: %e\n',var(logmodY));
    fprintf('DESVIO PADRAO...: %e\n',sqrt(var(logmodY)));
    fprintf('MAXIMO..........: %e\n',max(logmodY));
    fprintf('MINIMO..........: %e\n',min(logmodY));
    fprintf('#########################################################\n');
    fprintf('#########################################################\n');
    fprintf('ESTATISTICA DOs CAMPOs log-PERMEABILIDADE\n')
    fprintf('MEDIA...........: %e\n',mean(logK));
    fprintf('VARIANCIA.......: %e\n',var(logK));
    fprintf('DESVIO PADRAO...: %e\n',sqrt(var(logK)));
    fprintf('MAXIMO..........: %e\n',max(logK));
    fprintf('MINIMO..........: %e\n',min(logK));
    fprintf('#########################################################\n');
end
if(tipophi==1)
    fprintf('#############################\n')
    fprintf('Coeficientes para o campo phi\n')
    fprintf('  Valor de M: %f\n  Valor de S: %f\n',MEDIAPHI,STDPHI)
else
    fprintf('#############################\n')
    fprintf('Coeficientes para o campo phi\n')
    fprintf('  Valor de M: %f\n  Valor de S: %f\n',M,alpha)
end
fprintf('#############################\n')
fprintf('Coeficientes para o campo k\n')
fprintf('  Valor de M: %e\n  Valor de S: %f\n',mean(vetMM),mean(vetaalpha))
fprintf('#############################\n')
fprintf('Coeficientes para o campo E\n')
fprintf('  Valor de M: %e\n  Valor de S: %f\n',mean(veteM),mean(vetealpha))
fprintf('#############################\n')
fprintf('#############################\n')
fprintf('## MEDIA DE P (geral) #######\n')
fprintf('P_med = %f\n',MEDIAP/N)
fprintf('## DESVIO PADRAO DE P #######\n')
fprintf('SIGMA P_med = %f\n',sqrt(var(varTotalP)));
fprintf('## MAXIMO E MINIMO DE P #####\n')
fprintf('MAXIMO = %f\nMINIMO = %f\n',MAX,MIN);
fprintf('#############################\n')
fprintf('#############################\n')
fprintf('## MEDIA DE K (geral) #######\n')
fprintf('K_med = %e\n',MEDIAK/N)
fprintf('## DESVIO PADRAO DE K #######\n')
fprintf('SIGMA K_med = %e\n',sqrt(var(varTotalK)));
fprintf('#############################\n')
fprintf('#############################\n')
fprintf('## MEDIA DE E (geral) #######\n')
fprintf('E_med = %e\n',MEDIAE/N)
fprintf('## DESVIO PADRAO DE E #######\n')
fprintf('SIGMA E_med = %e\n',sqrt(var(varTotalE)));
fprintf('#############################\n')
sqrt(var(varTotalE))/(MEDIAE/N)
clear varTotalE varTotalK varTotalP vetMM vetaalpha veteM vetealpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(stat==1)
    a=max(porosidade);
    b=min(porosidade);
    x=[b*0.9:(a-b)/1000:a*1.1];
    y=0.0*x;
    yE=0.0*x;
    y =kcm(co,SS,x);
    if(self==1)
        for i=1:size(yE,2)
            yE(i) = selfcons(x(i),G,Ks,Kf,nu);
        end
    else
        yE= sc(SE,sY,x);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(8)
    figureKC(porosidade,permeabilidade,x,y,0)
    base=[home 'PoroXperm_' nome]
%    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    figureExphi(porosidade,young,x,yE,self,0)
    base=[home 'PoroXyoung_' nome]
%    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(9)
    if(tipophi==1)
        NORMAL(porosidade,mean(porosidade),sqrt(var(porosidade)),'Porosity')
    else
        LOGNORMAL(porosidade,mean(porosidade),sqrt(var(porosidade)),'Porosity')
    end
    base=[home 'PoroDist_' nome];
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(10)
    LOGNORMAL(permeabilidade,mean(permeabilidade),sqrt(var(permeabilidade)),'Permeability')
    base=[home 'PermDist_' nome];
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(11)
    LOGNORMAL(young,mean(young),sqrt(var(young)),'Youngs Modulus')
    base=[home 'YoungDist_' nome];
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(12)
    NORMAL(logK,mean(logK),sqrt(var(logK)),'Log-Permeability')
    base=[home 'LogPermDist_' nome];
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(14)
    NORMAL(logmodY,mean(logmodY),sqrt(var(logmodY)),'Log-Youngs Modulus')
    base=[home 'LogYoungDist_' nome];
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
end
cont;
clear
