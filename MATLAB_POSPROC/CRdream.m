clear;
d  = 5;
nCR= 3;
N  = 5;
delta = 2;
pg = 0.2;
cstar= 1e-12;
CR = [1:nCR]/nCR
c  = 0.1
pCR= ones(1,nCR)/nCR
[J,nid]=deal(zeros(1,nCR));
id = zeros(1,N);
%
x = nan(d,N);
X = randn(N,d);
x(1:d,1:N) = reshape(X',d,N);
lambda = unifrnd(-c,c,N,d);
dX = zeros(1,d);
for i=1:N
    a= 3;
    b= 2;
     D     = randsample([1:delta],1,'true')
     id(i) = randsample(1:nCR,1,'true',pCR)
     z     = rand(1,d)
     A     = find(z<CR(id(i)))
     d_star= numel(A)
     if(d_star == 0)
         [~,A] = min(z)
         d_star= 1
     end
     gamma = 2.38/sqrt(2*D*d_star)
     g     = randsample([gamma 1],1,'true',[1-pg pg])
     dX(A) = cstar*randn(1,d_star) + (1+lambda(i))*g*sum(X(a,A)-X(b,A),1)
 end
