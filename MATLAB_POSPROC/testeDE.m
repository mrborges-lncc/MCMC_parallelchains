clear
close all
d = 30;
NP= 5000;
X =randn(d,NP);
M = 1000;
x = X;
y = X;
beta = 1.0;
gamma = 2.38/sqrt(2*d);
alpha = 1.0/sqrt(2);
for i=1:M
%     x = randn(d,NP);
    for j=1:NP
        [a b]  = chooser(j,NP);
        y(:,j) = sqrt(beta)*(alpha)*x(:,j) + (alpha^2)*beta*(x(:,a)-x(:,b));
    end
    x = y;
end
mean(mean(y))
mean(var(y))