function [out] = load_data(file_ref,num_datatype,cut)
    for i = 1:num_datatype
        name   = strtrim(strjust(file_ref(i,:)));
        aux    = load(name);
        out{i} = aux(1+cut(1):end-cut(2),:);
    end
end

