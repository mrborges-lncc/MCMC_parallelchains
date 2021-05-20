function [C_kc, E0] = findKC_Spriggs_parameters(K,E,poro,spriggs)

%% Kozeny-Carman
syms f(x) x

phi  = poro;
f(x) = x * (phi^3)/(1-phi)^2 - K;

C_kc = double(solve(f,x))*0.85;

%phi = linspace(0.05,0.4,200);
%k   = C_kc * (phi.^3)./(1-phi).^2;

%plot(phi,k);

%% Spriggs relation

E0 = E / exp(-spriggs * poro);

%E = E0 * exp(-spriggs * phi);

%plot(phi,E)