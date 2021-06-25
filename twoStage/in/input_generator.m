%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='latex';
TOL=1e-7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {...
    'Experiment name:',...
    'Name of the samples:',...
    'Number of stages (1 or 2):',...
    'Number of priors $\theta$s:',...
    'Maximum number of simulation:',...
    'Maximum number of (different) selected fileds:',...
    'Data type number:'};
%
dlg_title = 'EXPERIMENT INPUT';
num_lines = [1 70];
def    = {'TwoPhase3D_DE','amostra','1','2','100000','5000','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
test   = int64(size(answer,1));
if(test==0)
    clear all
    error('ERROR');
end
ans = char(answer);
nome_exp    = char(answer(1,:));
nome_sample = char(answer(2,:));
ns   = int64(str2num(ans(3,:)));
np   = int64(str2num(ans(4,:)));
maxs = int64(str2num(ans(5,:)));
maxf = int64(str2num(ans(6,:)));
ndata= int64(str2num(ans(7,:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA MANAGEMENT
prompt = [];
nwells = [];
lkh  = [];
init = [];
final= [];
norma= [];
ef   = [];
ec   = [];
fref = [];
fsam = [];
for i=1:ndata
    dlg_title = ['Data set ' num2str(i)];
    prompt = {...
        'Number of wells:',...
        'Likelihood type:',...
        'Initial data:',...
        'Final data:',...
        'Normalize data? (yes or not):',...
        'Reference data file:',...
        'Sample data file:',...
        'Error at likelihood $\sigma_f$ (fine scale solution):'};
    def    = {'1','2','20','60','yes',...
        '../twophaseflow/exp/pres/pres_ref_0.dat',...
        '../twophaseflow/exp/pres/pres_amostra_0.dat',...
        '1.5e-05'};
    if ns == 2
        prompt = [prompt {'Error at likelihood $\sigma_c$ (coarse scale solution):'}];
        def    = [def {'2e-04'}];
    end
    num_lines = [1 70];
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test   = int64(size(answer,1));
    if(test==0)
        clear all
       error('ERROR');
    end
    ans = char(answer);
    nwells = [nwells; int64(str2num(ans(1,:)))];
    lkh    = [lkh; int64(str2num(ans(2,:)))];
    init   = [init; int64(str2num(ans(3,:)))];
    final  = [final; int64(str2num(ans(4,:)))];
    nor    = char(answer(5,:));
    nomer  = char(answer(6,:));
    nomes  = char(answer(7,:));
    fref   = char(fref,nomer);
    fsam   = char(fsam,nomes);
    ef     = [ef; double(str2num(ans(8,:)))];
    if ns == 2, 
        ec = [ec; double(str2num(ans(9,:)))];
    else
        ec = [ec; double(str2num(ans(8,:)))];
    end
    if nor == 'yes'
        norma = [norma; 1];
    else
        norma = [norma; 0];
    end
end
fref = fref(2:end,:);
fsam = fsam(2:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIORS
prop = [];
stcd = [];
njump= [];
fjump= [];
jump = [];
gera = [];
Tmat = [];
for i=1:np
    dlg_title = ['Prior ' num2str(i)];
    list   = {'gera_LABTRANGEORW3', 'FORTRAN', 'gera_LABTRANGEO',...
        'FORTRAN_RW', 'gera_LABTRANGEORW', 'gera_LABTRANGEORW2',...
        'MVN', 'MVN1', 'MVN2', 'FORTRAN_RW1', 'FORTRAN_RW2',...
        'FORTRAN_RW3D', 'FORTRAN_LBD','FORTRAN_KL3D','FORTRAN_KL3D_2'};
    [indx,tf] = listdlg('ListString',list);
    indx = indx-1;
    if indx == 9 , indx = 31; end
    if indx == 10, indx = 34; end
    if indx == 11, indx = 32; end
    if indx == 12, indx = 33; end
    gera = [gera; indx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt = {...
        'Proposal (2==RW; 3==AM; 4==DE 5==DREAM):',...
        'Stochastic dimension (KL):',...
        'Jump (0 == fixed; 1 == variable):',...
        'Jump frequency:',...
        'Jump size (mean or fixed):',...
        'Name of eigen-vector matrix:'};
    num_lines = [1 70];
    def    = {'4','13005','0','10','0.25','../gera_KL/MATLAB/out/avet1_510x510x20_51x51x5_l50x50x10_M13005.bin'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test   = int64(size(answer,1));
    if(test==0)
        clear all
       error('ERROR');
    end
    ans    = char(answer);
    prop   = [prop; int64(str2num(ans(1,:)))];
    stcd   = [stcd; int64(str2num(ans(2,:)))];
    njump  = [njump; int64(str2num(ans(3,:)))];
    fjump  = [fjump; int64(str2num(ans(4,:)))];
    jump   = [jump; double(str2num(ans(5,:)))];
    nome   = char(answer(6,:))
    Tmat   = char(Tmat,nome);
end
Tmat = Tmat(2:end,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list = {'Blackbox','SIMULADOR','simuladorRigido','simul_comp',...
    'UW SIMULATOR','SIMULADOR_VISCOELASTICO',...
    'twophaseflow','twophaseflowOCT'};
[indx,tf] = listdlg('ListString',list);
simul = int64(indx-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ns == 2
    prompt = [];
    for i=0:np-1
        name = ['Upscaling ' num2str(i) ':'];
        prompt = [prompt; {name}];
    end
    dlg_title = 'Upscaling';
    num_lines = [1 70];
    def    = {'0','1','2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test   = int64(size(answer,1));
    if(test==0)
        clear all
       error('ERROR');
    end
    ans = char(answer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    up = zeros(np,1);
    for i=1:np
        up(i) = int64(str2num(ans(i,:)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('entra_teste.in','w');
fprintf(fileID,'%10d\n',ns);
fprintf(fileID,'%10d\n',np);
for i=1:np
    fprintf(fileID,'%s\n',nome_sample);
end
if ns == 2
    for i=1:np
        fprintf(fileID,'%10d',up(i));
    end
else
    fprintf(fileID,'%10d',0);
end
fprintf(fileID,'\n');
fprintf(fileID,'%10d\n',maxs);
fprintf(fileID,'%10d\n',maxf);
fprintf(fileID,'%10d\n',simul);
fprintf(fileID,'%10d',ndata);
fprintf(fileID,'%10d\n',norma(1));
for i=1:ndata
    fprintf(fileID,'%10d%10d%10d%10d\n',nwells(i),lkh(i),init(i),final(i));
end
for i=1:np
    fprintf(fileID,'%10d%10d%10d%10d%10d%10d%12.5e\n',gera(i),prop(i),1,stcd(i),...
        njump(i),fjump(i),jump(i));
    fprintf(fileID,'%s\n',Tmat(i,:));
end
for i=1:ndata
    fprintf(fileID,'%12.5e\n',ef(i));
    fprintf(fileID,'%12.5e\n',ec(i));
end
for i=1:ndata
    fprintf(fileID,'%s\n',fref(i,:));
    fprintf(fileID,'%s\n',fsam(i,:));
end
fprintf(fileID,'%s\n',nome_exp);
NSEED = 33;
for i=1:NSEED
    fprintf(fileID,'%12d',int32(rand(1,1)*1e6));
end
fprintf(fileID,'\n%12d',0);
fclose(fileID);
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%