function [K, phi] = geologia(G,n,filenx,fileny,filenz,filephi,nD,beta,rho,...
    phibeta,phirho,depth,kc,ckc,heter,heterp,printa,nome)
    if heter == 1
        if kc
            Y   = load_perm(G,filephi,filephi,filephi,depth,n,nD);
            phi = phibeta .* exp(phirho * Y(:,1));
            clear Y
            %% Kozeny–Carman relation 
            perm = ckc * (phi.^3) ./ ((1.0 - phi).^2);
            sv   = 1e-03;
            perm = perm + lhsnorm(0,sv,G.cells.num).*perm;
            K    = [perm perm perm];
            clear perm
        else
            Y = load_perm(G,filenx,fileny,filenz,depth,n,nD);
            K = beta .* exp(rho * Y);
            phi = betaphi;
            clear Y
            if heterp == 1
                Yphi = load_perm(G,filephi,filephi,filephi,depth,n,nD);
                rhophi = 0.20;
                betaphi= phi;
                phi = betaphi .* exp(rhophi * Yphi(:,1));
                clear Yphi
            end
        end
        save('perm.mat','K');
        save('phi.mat','phi');
    end
    if heter == 2
        load('perm.mat','K');
        load('phi.mat','phi');
    end
    if heter ~= 2 && heter ~= 1
        if kc
            perm   = ckc * (phibeta.^3) ./ ((1.0 - phibeta).^2);
            K   = [perm perm perm];
            clear perm
        else
            K   = beta;
            phi = phibeta;
        end
    end
    fprintf('\nPorosity:\n Mean.....: %4.2f\n Std......: %4.2f\n Max......: %4.2f\n Min......: %4.2f\n',mean(phi),std(phi),max(phi),min(phi));
    fprintf('\nPermeability:\n Mean.....: %4.2e\n Std......: %4.2e\n Max......: %4.2e\n Min......: %4.2e\n',mean(K(:,1)),std(K(:,1)),max(K(:,1)),min(K(:,1)));
    if printa == 1 && heter == 1 && kc
        plotKC(phi,K(:,1),ckc)
        base=['figuras/perm_phi_KC_' nome]
        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r600', base);
        pause(3); clf; close all
    end
end

