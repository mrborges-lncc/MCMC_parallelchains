clear;
close all;
m = 20;
d = 1;
n = 100000;
b = n/20;
X=randn(n,d,m);
X=zeros(n,d,m);
for i=1:d
    X(:,i,:) = i+sqrt(i)*randn(n,m);
end
N=n/4;
for i=1:d
    X(1:N,i,:) = 100*rand(1,1)+sqrt(100*rand(1,1))*randn(N,m);
end
xcv = [];
cv1 = [];
W   = [];
V   = [];
for i=1:n
    if(mod(i,b)==0)
        ini = i/2;
        fim = i;
        xcv = [xcv; fim];
        [c,v,w] = convergeRMatrix(X(ini:fim,:,:));
        cv1 = [cv1; c];
        W   = [W; w];
        V   = [V; v];
    end
end
Rfigure(xcv,cv1)
Rfigure2(xcv,(W),(V))
