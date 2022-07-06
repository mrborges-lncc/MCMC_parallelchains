function [MMat, mutilde] = CONDmatrix(phi,lambda,M,n_dados,pnode,dados)
    corretor = 0;
    S=zeros(n_dados,n_dados);
    phicond=zeros(n_dados,M);
    for i = 1 : n_dados
        for j = 1 : n_dados
            for k = 1 : M
                S(i,j) = S(i,j) + ...
                    phi(pnode(i),k) * phi(pnode(j),k) * lambda(k)+corretor;
                phicond(i,k) = phi(pnode(i),k)+corretor;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R=zeros(M,n_dados);
    for i = 1 : M
        for k = 1 : n_dados
            R(i,k) = phi(pnode(k),i)+corretor;
        end
    end
    if(n_dados>0)
        LM=diag(lambda(1:M));
        S=R'*LM*(R);
    end
    %**************************************************************************
    if(n_dados>0)
        mutilde = ((LM).^(1/2))*R*inv(S)*dados;
        MMat=eye(M,M)-((LM).^(1/2))*R*inv(S)*R'*((LM).^(1/2));
    end
    return
end

