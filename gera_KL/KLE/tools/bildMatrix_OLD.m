function [C, sigma, num_elem] = bildMatrix_OLD(Lx,Ly,Lz,nx,ny,nz,...
    eta1,eta2,eta3,beta,nu,ntipo,varY)
    %** Parameters adjust *************************************************
    alpha    = 1.0;
    cutoff   = Lx/double(nx);
    sigma    = sqrt(varY);
    num_elem = nx*ny*nz;
    hx=Lx/double(nx);
    hy=Ly/double(ny);
    hz=Lz/double(nz);
    if(ntipo==2)
        ft = sqrt(2.0)*(varY);
        fat = (varY)*(hx^beta)*((hx/cutoff)^-beta);
    else
        fat=1.0;
        ft=varY;
    end
%
    coord=zeros(3,num_elem);
    k=0;
    for l=1:nz
        for j = 1:ny
            for i = 1:nx
                k = k + 1;
                coord(1,k) = (i-1)*hx+hx/2;
                coord(2,k) = (j-1)*hy+hy/2;
                coord(3,k) = (l-1)*hz+hz/2;
            end
        end
    end
    wP=1;
    %**************************************************************************
    % bulding the Covariance matrix *******************************************
    C = zeros(num_elem,num_elem);
    ephilon = hx;
    %
    if ntipo==2
        aux1 = (fat*alpha^(beta)*sigma^2*(hx/4)^2*wP');
        for ei = 1:num_elem
            zi = coord(:,ei)';
            C(ei,ei) =  alpha^(beta)*sigma^2*hx^2*sqrt(2);
            xv = [zi(1),zi(2),zi(3)];
            for ej = ei+1:num_elem
                zj = coord(:,ej)';
                yv = [zj(1),zj(2),zj(3)];
                aux = (ephilon./sqrt((xv(:,1)-yv(:,1)).^2+...
                    (xv(:,2)-yv(:,2)).^2+(xv(:,3)-yv(:,3)).^2)).^beta;
                C(ei,ej) = aux1*aux;
                C(ej,ei) = C(ei,ej);
            end
        end
    end
    if ntipo==1
        l1=eta1*eta1;
        l2=eta2*eta2;
        l3=eta3*eta3;
        for ei = 1:num_elem
            zi = coord(:,ei)';
            xv = [zi(1),zi(2),zi(3)];
            for ej = ei:num_elem
                zj = coord(:,ej)';
                yv = [zj(1),zj(2),zj(3)];
                 C(ei,ej) = Cov3D(xv,yv,varY,l1,l2,l3);
                 C(ej,ei) = C(ei,ej);
            end
        end
    end
    if ntipo==3
        l1=eta1^2;
        l2=eta2^2;
        l3=eta3^2;
        for ei = 1:num_elem
            zi = coord(:,ei)';
            xv = [zi(1),zi(2),zi(3)];
            for ej = ei:num_elem
                zj = coord(:,ej)';
                yv = [zj(1),zj(2),zj(3)];
                C(ei,ej) = Cov3Df(xv,yv,varY,l1,l2,l3);
                C(ej,ei) = C(ei,ej);
            end
        end
    end
    if ntipo==4
        l1=eta1*eta1;
        l2=eta2*eta2;
        l3=eta3*eta3;
        for ei = 1:num_elem
            zi = coord(:,ei)';
            xv = [zi(1),zi(2),zi(3)];
            for ej = ei:num_elem
                zj = coord(:,ej)';
                yv = [zj(1),zj(2),zj(3)];
                C(ei,ej) = matern(xv,yv,varY,l1,l2,l3,nu);
                C(ej,ei) = C(ei,ej);
            end
        end
    end
    clear xv yv wP aux coord zi zj

    return
end