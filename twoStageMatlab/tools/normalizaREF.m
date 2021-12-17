function [outmm, out] = normalizaREF(input,normaliza)
    sz = length(input);
    outmm = zeros(sz,2);
    if normaliza == 1
        for i = 1 : sz
            aux = input{i};
            outmm(i,:)   = [min(min(aux(:,2:end))) max(max(aux(:,2:end)))];
            aux(:,2:end) = (aux(:,2:end) - outmm(i,1))/...
                (outmm(i,2) - outmm(i,1));
            out{i} = aux;
        end
    else
        out = input;
        outmm(:,2) = 1.0;
    end
end

