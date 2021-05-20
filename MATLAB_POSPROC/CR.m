clear;
nCR= 3;
N  = 5;
CR = [1:nCR]
pCR= ones(1,nCR)/nCR
[J,nid]=deal(zeros(1,nCR))
id = zeros(1,N)
%
for i=1:N
  id(i) = randsample(1:nCR,1,'true',pCR)
end
