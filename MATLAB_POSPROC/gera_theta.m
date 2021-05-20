clear;
N=13005;
%A=randn(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_file = '../gera_KL/FORTRAN_KL3D/out/theta.dat';

for j=1:1
%    line_file = '../forecast/3Dfields/thetasMC/thetaMC_';
%    line_file = [line_file num2str(j-1) '.dat']
    fid = fopen(line_file,'w');
    A=lhsnorm(0.0,1.0,N);
    for i=1:N
        fprintf(fid,'%15.8e\n',A(i));
    end
    fclose(fid);
    MEDIA=mean(A)
    VAR=var(A)
end
    
clear