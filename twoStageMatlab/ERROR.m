clear all
close all
ini = 1;
fim = 10;
nt  = fim - ini + 1;
NC  = [1:5];
expname = 'DE';
homee   = './error/error';
er1     = [];
er2     = [];
%% load files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz  = 1e32;
for chain = NC
    name = [homee expname '_chain' num2str(chain,'%d') '.dat']
    erro = load(name);
    sz   = min(sz,size(erro,1))
end
for chain = NC
    name = [homee expname '_chain' num2str(chain,'%d') '.dat']
    erro = load(name);
    er1  = [er1 erro(1:sz,1)];
    er2  = [er2 erro(1:sz,2)];
end
plot(1:sz,er1)