function[Xi,Lx,Ly,nx,ny]= leitura2D(name)
%
        fid = fopen(name,'r');
        mattamp = fscanf(fid,'%f');
        fclose(fid);
        inf = mattamp(1:4);
        Lx = inf(1);
        Ly = inf(2);
        nx = inf(3);
        ny = inf(4);
        dx = Lx/nx;
        dy = Ly/ny;
        Xi=zeros(nx*ny,1);
    %
        mattamp = mattamp(9:length(mattamp));
%         permmap=zeros(ny,nx);
        k=0;
        loc=0;
        for j=ny:-1:1
            k=k+1;
            if(mattamp(k)~=ny-j)
                disp('erro1')
                break
            end
            for i=1:nx
                k=k+1;
%                 permmap(j,i)=mattamp(k);
                loc = loc + 1;
                Xi(loc) = mattamp(k);
            end
            k=k+1;
            if(mattamp(k)~=192837465)
                disp('erro2')
                break
            end       
        end
        clear mattamp
end
