function imprime(Lx,Ly,nx,ny,ntipo,beta,Xi,nr,home,name,iprt)
    nfile= num2str(nr,5);
    if iprt==2
        bb2 = [home name '_' nfile '.mat']    
        save(bb2,'Xi');
    else
       for dim = 1:2
           bb2 = [home name '_d' num2str(dim) '_' nfile '.dat']
           outfile2 = fopen(bb2, 'wt');
           if(iprt==1)
                fprintf(outfile2,'%f\n',Lx);
                fprintf(outfile2,'%f\n',Ly);
                fprintf(outfile2,'%d\n',nx);
                fprintf(outfile2,'%d\n',ny);
                fprintf(outfile2,'%d\n',ntipo);
                fprintf(outfile2,'%f\n',beta);
                fprintf(outfile2,'%d\n',2);
                fprintf(outfile2,'%d\n',2);
                el=0;
                for i=1:ny
                  fprintf(outfile2,'%d\n',i-1);
                  for j=1:nx
                      el=el+1;
                      fprintf(outfile2,'%f ',Xi(el,dim));
                  end
                  fprintf(outfile2,'\n192837465\n');
                end
           end
           if(iprt==0)
               fprintf(outfile2,'# campos KL\n');
               fprintf(outfile2,'%d x %d x 1\n',nx,ny);
               for i=1:nx*ny
                   fprintf(outfile2,'%f\n',Xi(i,dim));
               end
           end
           if(iprt==3)
               aux = zeros(ny,nx);
               for i=1:ny
                   ini = (i-1)*nx+1;
                   aux(i,:) = Xi(ini:ini+nx-1,dim);
               end
           end
           fclose(outfile2);
       end
    end
end
