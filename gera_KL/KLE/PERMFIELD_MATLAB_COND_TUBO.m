% Single phase flow Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 25/03/2020
% Based on MRST 2019
% Author: Marcio Borges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

tstart =  tic;
addpath ./tools/
addpath ~/Dropbox/mrst-2022a/

startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 36.0;
Ly  = 36.0;
Lz  = 50.0;
nx  = 36;
ny  = 36;
nz  = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printafig = 1;
printaMean= 1;
tipo_prt  = 0;
depth     = 0e3;
condit    = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ini = 0;
fim = 3;
NT  = fim - ini + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 0.817651;
beta= 6.6707e-14;      %% Factor to permeability
nome= 'permREF';
variav = '\kappa';
%
rho = 0.125;
beta= 0.275;
nome= 'phiREF';
variav = '\phi';
%
% rho    = 0.457;
% beta   = 1.0225e10;
% nome   = 'EREF';
% variav = '\mathsf{E}';
variav = '';
nome  = 'TUBsexpCOND'
homef = './figuras/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FILENAME, PATHNAME] =uigetfile({'~/Dropbox/KLE/fields/se*.dat'}, 'LOAD DATA');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
filen = filen(1:k);
nini = int32(str2num(filem(k+1:end)));
verbose = true;
nD = '3D';
color = 'none';
%color = 'k';
vw  = [35 20];
%vw  = [0 90];
et = 3;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if condit
    dados = load('~/Dropbox/KLE/in/input_cond.dat');
    dados(:,3) = Lz - dados(:,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx  = Lx/double(nx);
dy  = Ly/double(ny);
dz  = Lz/double(nz);
G   = cartGrid([nx ny nz],[Lx Ly Lz]*meter^3);
G.nodes.coords(:, 3) = Lz - G.nodes.coords(:, 3);
% G.nodes.coords(:, 2) = G.nodes.coords(:, 2);
% G.nodes.coords(:, 1) = G.nodes.coords(:, 1);
G   = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx2 = dx/2; dy2 = dy/2; dz2 = dz/2;
if condit
    cells = [];
    for j = 1 : size(dados,1)
        v = dados(j,1:3);
        for i = 1 : G.cells.num
            c = G.cells.centroids(i,:);
            if c(1)-dx2 <= v(1) && v(1) <= c(1)+dx2
                if c(2)-dy2 <= v(2) && v(2) <= c(2)+dy2
                    if c(3)-dz2 <= v(3) && v(3) <= c(3)+dz2
                        cells = [cells; i];
                    end
                end
            end
        end
    end
    Ge = cell(size(dados,1),1);
    for j = 1 : size(dados,1)
        c  = G.cells.centroids(cells(j),:);
        Gb = cartGrid([1 1 1],[dx dy dz]);
        Gb.nodes.coords(:, 3) = c(3) + Gb.nodes.coords(:, 3)-dz2;
        Gb.nodes.coords(:, 2) = c(2) + Gb.nodes.coords(:, 2)-dy2;
        Gb.nodes.coords(:, 1) = c(1) + Gb.nodes.coords(:, 1)-dx2;
        Ge{j} = computeGeometry(Gb);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center = [Lx/2 Ly/2 Lz/2];
raio   = 18.0;
r = zeros(nx * ny * nz,1);
c = G.cells.centroids;
for i = 1 : nx * ny * nz
    r(i) = norm(c(i,1:2) - center(1:2));
end
GC = removeCells(G, r > raio);
p  = find(r <= raio);
% plotGrid(GC); axis equal tight; view(vw); box 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOLOGIC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanYfield = zeros(G.cells.num,1);
for n = ini:fim
    snum = num2str(n,'%d');
    Yorig = load_perm(G,filen,filen,filen,depth,n,nD,tipo_prt,nx,ny,nz);
    Yorig = Yorig(:,1);
    Y = Yorig(p,1);
    K = beta .* exp(rho * Y);
    meanYfield = meanYfield + Yorig;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rock = makeRock(GC, K, 0.20);
    mK   = mean(rock.perm/(milli*darcy));
    sK   = std(rock.perm/(milli*darcy));
    fprintf('\n==============================================================\n')
    fprintf('Mean Y....: %4.1f    \t | \t std Y....: %4.1f   \n',mean(Y),std(Y));
    fprintf('Mean K....: %4.1f mD \t | \t std K....: %4.1f mD\n',mK,sK);
    fprintf('==============================================================\n')
    clear K
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if printafig == 1
        lim = [-3. 3.];
        %plot_rock(Y,G,'Yn','$Y$',color,lim,vw,1);
        plot_rock_poroZinv(rock.perm(:,1),GC,'Y',beta,rho,...
            ['$\mathsf{Y}_{' variav '}\!\left( \mathbf{x}  \right)$'],color,lim,vw,20);
        if condit
            hold on
            for j = 1 : size(dados,1)
                plotGrid(Ge{j},'FaceColor','none','EdgeColor','k',...
                    'LineWidth',2);
                fprintf('\n==============================================')
                fprintf('\nPoint (%3.2f, %3.2f, %3.2f) = %3.2f',...
                    dados(j,1:3),Yorig(cells(j)))
                fprintf('\n==============================================')
            end
        end
        base=[homef 'Y_' nome '_' snum];
        set(gcf,'PaperPositionMode','auto');
    %     print('-depsc','-r300', base);
        print('-dpng' ,'-r300', base);
        pause(2); close all
    % %     lim = [(mean(rock.perm(:,1)) - 0.0125*std(rock.perm(:,1)))...
    % %         (mean(rock.perm(:,1)) + 0.25*std(rock.perm(:,1)))];
    % %     plot_rock(rock.perm(:,1),G,'Yn',variav,color,lim,vw,2);
    %     plot_rock_poro(rock.perm(:,1),G,'Yn',beta,rho,['$' variav '$'],...
    %         color,lim,vw,12);
    %     base=[homef nome '_' snum];
    %     set(gcf,'PaperPositionMode','auto');
    %     print('-depsc','-r600', base);
    %     pause(1); close all
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printaMean == 1
    meanYfield = meanYfield / NT;
    lim = [min(meanYfield) max(meanYfield)];
    plot_rock_poroZinv(meanYfield(p,1),GC,'Yn',beta,rho,...
        ['$\langle\mathsf{Y}_{' variav '}\!\left( \mathbf{x}  \right)\rangle$'],color,lim,vw,30);
    if condit
        hold on
        for j = 1 : size(dados,1)
            plotGrid(Ge{j},'FaceColor','none','EdgeColor','k',...
                'LineWidth',2);
            fprintf('\n==============================================')
            fprintf('\nPoint (%3.2f, %3.2f, %3.2f) = %3.2f',...
                dados(j,1:3),meanYfield(cells(j),1))
            fprintf('\n==============================================')
        end
    end
    base=[homef 'meanY_' nome '_' snum];
    set(gcf,'PaperPositionMode','auto');
    %     print('-depsc','-r300', base);
    print('-dpng' ,'-r300', base);
    %     pause(1); close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
