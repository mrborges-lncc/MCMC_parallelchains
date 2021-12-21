function T = load_KL(file_KL,num_rockpar,fine_mesh,KLM)  
    dim = fine_mesh(1) * fine_mesh(2) * fine_mesh(3);
    for i = 1:num_rockpar
        name   = strtrim(strjust(file_KL(i,:)));
        fileID = fopen(name,'r');
        TKL    = fread(fileID,'single');
        M      = int64(size(TKL,1)/dim);
        if KLM > M
            error('Error at stochastic dimention')
        end
        TKL    = reshape(TKL,[dim,M]);
        T{i}   = TKL(:,1:KLM);
        fclose(fileID);
    end 
end

