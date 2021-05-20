function imprime3D(Lx,Ly,Lz,nx,ny,nz,ntipo,beta,Xi,nr,home,name,iprt)
   nfile= num2str(nr,'%d');
   if(iprt==2)
        bb2 = [home name '_' nfile '.mat']    
        save(bb2,'Xi');
   else
        for dim = 1:3
           if(iprt==1)
               bb2 = [home name '_d' num2str(dim) '_' nfile '.dat']
               outfile2 = fopen(bb2, 'wt');      
                fprintf(outfile2,'%f\n',Lx);
                fprintf(outfile2,'%f\n',Ly);
                fprintf(outfile2,'%f\n',Lz);
                fprintf(outfile2,'%d\n',nx);
                fprintf(outfile2,'%d\n',ny);
                fprintf(outfile2,'%d\n',nz);
                fprintf(outfile2,'%d\n',ntipo);
                fprintf(outfile2,'%f\n',beta);
                fprintf(outfile2,'%d\n',2);
                fprintf(outfile2,'%d\n',2);

                m=0;
                for z=1:nz
                    fprintf(outfile2,'%d\n',z-1);
                    for i=1:ny
                        fprintf(outfile2,'%d\n',i-1);
                        for j=1:nx
                            m=m+1;
                            fprintf(outfile2,'%g ',Xi(m,dim));
                        end
                        fprintf(outfile2,'\n192837465\n');
                    end
                end
                fclose(outfile2);
           end
           if(iprt==0)  
               bb2 = [home name '_d' num2str(dim) '_' nfile '.dat']
               outfile2 = fopen(bb2, 'wt');  
               fprintf(outfile2,'# campos KL\n');
               fprintf(outfile2,'%d x %d x %d\n',nx,ny,nz);
               for i=1:nx*ny*nz
                   fprintf(outfile2,'%12.5g\n',Xi(i,dim));
               end
               fclose(outfile2);
           end
        end
   end
end
%
