function [out] = normaliza(input,outmm)
    sz = length(input);
    for i = 1 : sz
        aux = input{i};
        aux(:,2:end) = (aux(:,2:end) - outmm(i,1))/(outmm(i,2) - outmm(i,1));
        out{i} = aux;
    end
    out = out';
end

