function [out] = load_data(file_ref,num_datatype)
    for i = 1:num_datatype
        name   = strtrim(strjust(file_ref(i,:)));
        out{i} = load(name);
    end
end

