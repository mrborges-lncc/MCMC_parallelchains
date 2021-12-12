clear all
close all
addpath ./tools/

ini = 1;
fim = 10;
nt  = fim - ini + 1;
NC  = [1:2];
expname = 'CW';
stage   = 1;
d       = 10000;
home    = '~/twoStageMatlab/';
home    = './';
homee   = [home 'error/error'];
homer   = [home 'out/restart'];
er1     = [];
er2     = [];
%% load files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz  = 1e32;
for chain = NC
    name = [homee expname '_chain' num2str(chain,'%d') '.dat'];
    erro = load(name);
    sz   = min(sz,size(erro,1));
end
for chain = NC
    name = [homee expname '_chain' num2str(chain,'%d') '.dat'];
    erro = load(name);
    er1  = [er1 erro(1:sz,1)];
    er2  = [er2 erro(1:sz,2)];
end

errorfigure(1:sz,er1,'$||\mathsf{F}_{\!\mathsf{ref}_{\ \! 1}} - \mathsf{F}_{\!\mathsf{s}_{\ \! 1}}||^{2}$')
name = ['figuras/ErrorF1_' expname];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600',name);

errorfigure(1:sz,er2,'$||\mathsf{F}_{\!\mathsf{ref}_{\ \! 2}} - \mathsf{F}_{\!\mathsf{s}_{\ \! 2}}||^{2}$')
name = ['figuras/ErrorF2_' expname];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600',name);

TOL = 1.0e-08;
counter = ones(length(NC),1);
fprintf('\n========================================\n')
for k = 1 : length(NC)
    chain= NC(k);
    name = [homee expname '_chain' num2str(chain,'%d') '.dat'];
    erro = load(name);
    aux  = erro(1,1);
    for i = 2 : size(erro,1)
        if abs(aux - erro(i,1)) > TOL
            aux = erro(i,1);
            counter(k,1) = counter(k,1) + 1;
        end
    end
    fprintf('Acceptance rate in chain %d.....: %4.2f\n',chain,100*counter(k,1)/size(erro,1))
end
fprintf('========================================\n')
chain = length(NC);
fprintf('\n========================================\n')
if stage == 2
    nome = [homer expname '_stage' num2str(stage,'%d') '_NC' ...
        num2str(chain,'%d') '_d' num2str(d,'%d') '_ccounter.mat'];
    ccounter = load(nome,'-mat');
    for k = 1 : length(ccounter.ccounter)
        fprintf('Acceptance rate in coarse scale problem chain %d.....: %4.2f\n',NC(k),100*ccounter.ccounter(k)/size(erro,1));
    end
end    
if stage == 2
    nome = [homer expname '_stage' num2str(stage,'%d') '_NC' ...
        num2str(chain,'%d') '_d' num2str(d,'%d') '_counter.mat'];
    counter = load(nome,'-mat');
    for k = 1 : length(counter.counter)
        fprintf('Acceptance rate in fine scale problem chain %d.....: %4.2f\n',NC(k),100*counter.counter(k)/size(erro,1));
    end
end
clear