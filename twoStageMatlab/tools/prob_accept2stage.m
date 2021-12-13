function [out] = prob_accept2stage(dataref,samplen,sample,sigma,...
    csamplen,csample,sigmac,num_datatype,cpostratio)
    aux = 0.0;
    auxn= 0.0;
    for k = 1 : num_datatype
        ref   = dataref{k};
        samp  = sample{k};
        sampn = samplen{k};
        csamp = csample{k};
        csampn= csamplen{k};
        aux    = aux - ((norm(ref(:,2:end) - samp(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        auxn   = auxn - ((norm(ref(:,2:end) - sampn(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        caux   = caux - ((norm(ref(:,2:end) - csamp(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        cauxn  = cauxn - ((norm(ref(:,2:end) - csampn(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
        
    end
    out = min(1,exp(aux - auxn + cauxn - caux));
    return    
end

