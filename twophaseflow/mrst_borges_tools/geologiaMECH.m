function [K, phi, E, nu, alpha, Bulkm, CompR] = geologiaMECH(G,n,filenx,...
    fileny,filenz,nD,K,vK,beta,rho,phibeta,phirho,E,spriggs,nu,alpha,...
    depth,kc,heter,printa,nome)
    if heter == 1
        if kc
            [ckc, E0] = findKC_Spriggs_parameters(K(1),E,phibeta,spriggs);
            Y   = load_perm(G,filenx,fileny,filenz,depth,n,nD);
            Y   = Y - mean(Y);
            Y   = (1.0./std(Y)).*Y;
            phi = phibeta * exp(phirho * Y(:,1));
            clear Y
            %% Kozeny–Carman relation 
            perm = ckc * (phi.^3) ./ ((1.0 - phi).^2);
            sv   = 1e-03;
            perm = perm + lhsnorm(0,sv,G.cells.num).*perm;
            K    = [perm perm perm];
            clear perm
            %% Sprigg’s representation
            E  = E0 * exp(-spriggs * phi);
            sv = 3e-05;
            E  = E + lhsnorm(0,sv,G.cells.num).*E;
            %% Bulk modulus
            Bulkm = (E)./(3.0*(1.0 - 2*nu));
            %% poisson ratio
            nu    = repmat(nu, G.cells.num, 1);
            %% Biot 
            alpha = repmat(alpha, G.cells.num, 1);
            %% Rock compressibility
            CompR = 1./Bulkm;
        else
            Y     = load_perm(G,filenx,fileny,filenz,depth,n,nD);
            Y     = Y - mean(Y);
            Y     = (1.0./std(Y)).*Y;
%             muY   = log(K) - 0.5 * log(vK./K.^2 + 1);
%             varY  = log(vK./K.^2 + 1);
%             K     = exp(sqrt(varY).*Y - muY)
            K     = beta * exp(rho * Y);
            Bulkm = E./(3.0*(1.0 - 2.0*nu));        %% Bulk modulus
            E     = repmat(E, G.cells.num, 1);
            nu    = repmat(nu, G.cells.num, 1);
            alpha = repmat(alpha, G.cells.num, 1);
            CompR = 1.0/Bulkm;                      %% Rock compressibility
            phi   = phibeta;
        end
    else
        if kc
            [ckc, E0] = findKC_Spriggs_parameters(K(1),E,phibeta,spriggs);
            phi  = repmat(phibeta, G.cells.num, 1);
            %% Kozeny–Carman relation 
            perm   = ckc * (phi.^3) ./ ((1.0 - phi).^2);
            K   = [perm perm perm];
            clear perm
            %% Sprigg’s representation
            E   = E0 * exp(-spriggs * phi);
            %% Bulk modulus
            Bulkm = (E)./(3.0*(1.0 - 2*nu));
            %% poisson ratio
            nu    = repmat(nu, G.cells.num, 1);
            %% Biot 
            alpha = repmat(alpha, G.cells.num, 1);
            %% Rock compressibility
            CompR = 1./Bulkm;
        else
            if heter == 2
                load('out/perm.mat','K');
            else
                K = K;
            end
            E     = repmat(E, G.cells.num, 1);
            nu    = repmat(nu, G.cells.num, 1);
            alpha = repmat(alpha, G.cells.num, 1);
            Bulkm = E./(3.0*(1.0 - 2.0*nu));        %% Bulk modulus
            CompR = 1.0/Bulkm;                      %% Rock compressibility
            phi   = phibeta;
        end
    end
    fprintf('\nPorosity:\n Mean.....: %4.2f\n Std......: %4.2f\n Max......: %4.2f\n Min......: %4.2f\n',mean(phi),std(phi),max(phi),min(phi));
    fprintf('\nPermeability:\n Mean.....: %4.2e\n Std......: %4.2e\n Max......: %4.2e\n Min......: %4.2e\n',mean(K(:,1)),std(K(:,1)),max(K(:,1)),min(K(:,1)));
    fprintf('\nYoung modulus:\n Mean.....: %4.2e\n Std......: %4.2e\n Max......: %4.2e\n Min......: %4.2e\n',mean(E),std(E),max(E),min(E));
    fprintf('\nBulk modulus:\n Mean.....: %4.2e\n Std......: %4.2e\n Max......: %4.2e\n Min......: %4.2e\n',mean(Bulkm),std(Bulkm),max(Bulkm),min(Bulkm));
    fprintf('\nRock compressibility:\n Mean.....: %4.2e\n Std......: %4.2e\n Max......: %4.2e\n Min......: %4.2e\n',mean(CompR),std(CompR),max(CompR),min(CompR));
    if printa == 1 && heter == 1 && kc
        plotKC(phi,K(:,1),ckc)
        base=['figuras/perm_phi_KC_' nome]
        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r600', base);
        pause(3); clf; close all
        plotEphi(phi,E,E0)
        base=['figuras/E_phi_spriggs_' nome]
%        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r600', base);
        pause(3); clf; close all
    end
end

