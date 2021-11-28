function [out] = erromediorel(cref,csample,num_datatype)
    out = 0.0;
    for k = 1 : num_datatype
        ref    = cref{k};
        sample = csample{k};
        out    = out + ((norm(ref(:,2:end) - sample(:,2:end)))/...
            (norm(ref(:,2:end))));
    end
    out = out^2;
    return
end