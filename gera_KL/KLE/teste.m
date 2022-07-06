clear all
M = 20000;
num_elem = M;
theta = single(lhsnorm(0,1,M));
Xi    = single(zeros(num_elem,1));
mY    = single(zeros(num_elem,1));
phi   = single(rand(M,M));
tSgera = tic;
    for el=1:num_elem
        Xi(el) =  mY(el) + phi(el,:)*theta;
    end
tEgera=toc(tSgera);
disp(['Tempo total gasto na geraccao dos campos: ' num2str(tEgera) ' seg.'])
disp('------------------------------');

tSgera = tic;
X = mY + sum(phi.*theta',2);
tEgera=toc(tSgera);
disp(['Tempo total gasto na geraccao dos campos: ' num2str(tEgera) ' seg.'])
disp('------------------------------');

[Xi(1:10) X(1:10); Xi(end-10:end) X(end-10:end)]
