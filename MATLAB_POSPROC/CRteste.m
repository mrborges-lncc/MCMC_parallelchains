clear;
close all
nc = 10;
cr = [];
p  = [];
for i=1:nc
    cr = [cr; double(i)/double(nc)]
    p  = [p; 1.0/double(nc)]
end
id=[];
M = 50000;
for i=1:M
%    id = [id; multinom(nc,p)];
  id = [id;unifdisc(nc,1)];
end
cont = zeros(nc,1);
for i=1:length(id)
    for j=1:nc
        if(id(i)==j)
            cont(j)=cont(j)+1;
            break
        end
    end
end
cont/M
hist(id)