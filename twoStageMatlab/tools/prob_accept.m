function [out] = prob_accept(dataref,samplen,sample,sigma,...
    num_datatype,coarse_post_ratio)
    aux = 0.0;
    auxn= 0.0;
    for k = 1 : num_datatype
        ref    = dataref{k};
        sample = sample{k};
        samplen= samplen{k};
        aux    = aux - ((norm(ref(:,2:end) - sample(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        auxn   = auxn - ((norm(ref(:,2:end) - samplen(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        
    end
    out = min(1,exp(aux - auxn));
    return    
end

