clear all; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NC = 3;
d  = 1;
Ni = [1];
Nf = [600];
M  = 20;
Nt = Nf - Ni + 1;
expname = 'RW';
home    = '~/twoStageMatlab/';
home    = './';
homet   = [home 'thetas/theta'];
vari    = '1';
X = zeros(Nt,d,NC);
%% READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for chain = 1 : NC
    k = 0;
    for n = Ni : Nf
        name = [homet expname '_chain' num2str(chain,'%d') '_t' vari...
            '_' num2str(n,'%d') '.dat']
        aux = load(name,'-mat');
        k = k + 1;
        X(k,:,chain) = aux.t(1:d);
    end
end

% for chain = 1 : NC
%     X(:,:,chain) = mvnrnd(zeros(d,1),eye(d),Nf);
% end
%X(:,:,NC) = 0.5 + 3*X(:,:,NC);

x = [];
R = [];
for i = M : M : Nt
    [r,v,w] = convergeRMatrix(X(i/2:i,:,:));
    x = [x; i];
    R = [R; r];
end
% Rfigure(x,R,R,R)
plot(x,R)