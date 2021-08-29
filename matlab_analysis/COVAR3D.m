clear all; close all; clc;
addpath ~/Dropbox/mrst_borges_tools/
addpath ~/Dropbox/mrst-2021a/
startup

color = 'none';
lim = [0 0];
vw  = [-35 20];
home = '~/fields/campos/';
home_fig = '../figuras/';
nome     = 'E';
ini = 0;
fim = 1999;
N   = fim - ini + 1;
cutr= 0.4; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 510.0;
Ly  = 510.0;
Lz  = 20.0;
nx  = 51;
ny  = 51;
nz  = 5;
depth = 1e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx  = Lx/double(nx);
dy  = Ly/double(ny);
dz  = Lz/double(nz);
G   = cartGrid([nx ny nz],[Lx Ly Lz]*meter^3);
G.nodes.coords(:, 3) = depth + G.nodes.coords(:, 3)*meter;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*meter;
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*meter;
G   = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filen = [home nome '_'];
cx = floor(nx/2);
cy = floor(ny/2);
cz = floor(nz/2);
covx = zeros(cx+1,1);
covy = zeros(cy+1,1);
covz = zeros(cz+1,1);
conx = zeros(cx+1,1);
cony = zeros(cy+1,1);
conz = zeros(cz+1,1);
disx = dx*[0:size(covx)-1].';
disy = dy*[0:size(covy)-1].';
disz = dz*[0:size(covz)-1].';

Z = zeros(N,nx,ny,nz);
k = 0;
for i = ini : fim
    n = num2str(i,'%d');
    Y = load_perm(G,filen,filen,filen,depth,i,'3D');
    Y = Y(:,1);
    k = k + 1;
    Z(k,:,:,:) = reshape(Y,nx,ny,nz);
end
meanY = mean(Z,1);
stdY  = std(Z,0,[1]);
varY  = var(Z,0,[1]);

zk = 0;
for i = ini : fim
    zk = zk + 1;
    fprintf('\ndone %d/%d',zk,N);
    %% cov x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = 0;
    for n = 0:cx
        m = m + 1;
        for k = 1:nz
            for j = 1:ny
                for d = 1:nx-n
                    conx(m) = conx(m) + 1;
                    covx(m) = covx(m) + (Z(zk,d,j,k) - meanY(1,d,j,k)) * ...
                    (Z(zk,d+n,j,k) - meanY(1,d+n,j,k));
                end
            end
        end
    end
    %% cov y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = 0;
    for n = 0:cy
        m = m + 1;
        for k = 1:nz
            for j = 1:nx
                for d = 1:ny-n
                    cony(m) = cony(m) + 1;
                    covy(m) = covy(m) + (Z(zk,j,d,k) - meanY(1,j,d,k)) * ...
                    (Z(zk,j,d+n,k) - meanY(1,j,d+n,k));
                end
            end
        end
    end
    %% cov z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = 0;
    for n = 0:cz
        m = m + 1;
        for k = 1:ny
            for j = 1:nx
                for d = 1:nz-n
                    conz(m) = conz(m) + 1;
                    covz(m) = covz(m) + (Z(zk,j,k,d) - meanY(1,j,k,d)) * ...
                    (Z(zk,j,k,d+n) - meanY(1,j,k,d+n));
                end
            end
        end
    end
%     plot_rock_poro(Y,G,'Yn',1,1,['$Y$'],...
%         color,lim,vw,2);
end
covx = covx./conx;
covy = covy./cony;
covz = covz./conz;
%% figure x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = int64(cutr * cx);
x  = disx(1:sz);
y  = log(covx(1:sz));
aux= corrcoef([x y]);
coefr=aux(1,2)*aux(1,2);
RL = polyfit(x,y,1);
CR=num2str(coefr,'%1.2e')
aa=num2str(-1/RL(1,1),'%1.2f')
bb=num2str(exp(RL(1,2)),'%1.2f');
xx = [0:(max(disx)-min(disx))/300:max(disx)].';
pp = exp(RL(1,1)*xx + RL(1,2));
yname = '$\mathcal{C}\left( {x} \right)$';
xname = '${x}$'
dat1  = 'measured';
dat2  = ['$\hat{\mathcal{C}}\left( {x} \right) =' bb '\exp\left(-\frac{x}{' aa '}\right)$']

covplot(disx, xx, covx, pp, xname, yname, dat1, dat2, 43)
base=[home_fig 'Covx_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);

%% figure y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = int64(cutr * cy);
x  = disy(1:sz);
y  = log(covy(1:sz));
aux= corrcoef([x y]);
coefr=aux(1,2)*aux(1,2);
RL = polyfit(x,y,1);
CR=num2str(coefr,'%1.2e')
aa=num2str(-1/RL(1,1),'%1.2f')
bb=num2str(exp(RL(1,2)),'%1.2f');
xx = [0:(max(disy)-min(disy))/300:max(disy)].';
pp = exp(RL(1,1)*xx + RL(1,2));
yname = '$\mathcal{C}\left( {y} \right)$';
xname = '${y}$'
dat1  = 'measured';
dat2  = ['$\hat{\mathcal{C}}\left( {y} \right) =' bb '\exp\left(-\frac{y}{' aa '}\right)$']

covplot(disy, xx, covy, pp, xname, yname, dat1, dat2, 44)
base=[home_fig 'Covy_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);

%% figure y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutr = 1.0;
sz = int64(cutr * cz);
x  = disz(1:sz);
y  = log(covz(1:sz));
aux= corrcoef([x y]);
coefr=aux(1,2)*aux(1,2);
RL = polyfit(x,y,1);
CR=num2str(coefr,'%1.2e')
aa=num2str(-1/RL(1,1),'%1.2f')
bb=num2str(exp(RL(1,2)),'%1.2f');
xx = [0:(max(disz)-min(disz))/300:max(disz)].';
pp = exp(RL(1,1)*xx + RL(1,2));
yname = '$\mathcal{C}\left( {z} \right)$';
xname = '${z}$'
dat1  = 'measured';
dat2  = ['$\hat{\mathcal{C}}\left( {z} \right) =' bb '\exp\left(-\frac{z}{' aa '}\right)$']

covplot(disz, xx, covz, pp, xname, yname, dat1, dat2, 45)
base=[home_fig 'Covz_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
