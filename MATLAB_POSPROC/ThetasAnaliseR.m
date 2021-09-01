clear;
close all
jump=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nthetas = 2;
M       = 1;
Nch_ini = 0;
Nch_fim = 5;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(0, 1, Nchains);
Nfim = [18 20 18 21 13 20];
%Nfim = repmat(900, 1, Nchains);
Nfim = Nfim(Nch_ini+1:Nch_fim+1)-2;
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
nome = 'TwoPhase3D_onlyPerm_RW_RK';
nome = 'TwoPhase3D_RW_RK';
base_name = ['prod_D2_' nome];
hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
% hom = '~/Dropbox/PROJETO_MCMC_RIGID/MCMCrw_onlyPerm/';
% hom = '../';
homf= '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
home = [hom 'twoStage/select_thetas/'];
%
filename = [];
%
for i=1:Nthetas
    filename = [filename;...
        [home 'theta_v' num2str(i,'%d') '_' nome]];
end
%
total=0;
pmax = -1e10;
pmin = 1e10;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    pchains = load([home '../out/nchain_' nome n '.dat']);
    sz = size(pchains,1) + 1;
    total = total + sum(pchains(Nini(j)+1:Nfim(j)+1,2));
    pmax  = max(pmax,sum(pchains(Nini(j)+1:Nfim(j)+1,2)));
    pmin  = min(pmin,sum(pchains(Nini(j)+1:Nfim(j)+1,2)));
end
total = int64(total);
pmin  = int64(pmin);
data1 = zeros(pmin,M,Nchains);
if(Nthetas == 2), data2 = zeros(pmin,M,Nchains); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:Nthetas
    file_name = [filename(k,:)];
    for j=1:Nchains
        n = num2str(chains(j),'%d');
        name = [home '../out/nchain_' nome n '.dat'];
        pchains = load(name);
        cont = 0;
        for i=Nini(j):Nfim(j)
            istr=num2str(i,5);
            fname = [file_name n '_' istr '.dat'];
            dat = load(fname);
            m   = pchains(i+1,2);
            for nc = 1:m
                cont = cont + 1;
                if cont > pmin, break; end
                if(k == 1)
                    data1(cont,:,j) = lhsnorm(0,1,M).';%dat(1:M);
                else
                    data2(cont,:,j) = dat(1:M);
                end
            end
        end
    end
end
%
[cv1,v1,w1] = convergeRMatrix(data1)
[cv2,v2,w2] = convergeRMatrix(data2)
fprintf('\n=========================================')
fprintf('\n=========================================')
fprintf('\nNumero total de dados: %d | %d\n',cont,total);
sz = size(data1);
x  = zeros(sz(1)*sz(2)*sz(3),1);
y  = zeros(sz(1)*sz(2)*sz(3),1);
k  = 0;
for i=1:sz(3)
    for j=1:sz(1)
        x(k*M+1:(k+1)*M) = data1(j,:,i).';
        y(k*M+1:(k+1)*M) = data2(j,:,i).';
        k = k + 1;
    end
end

