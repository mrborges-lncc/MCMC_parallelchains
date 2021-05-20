clear;
base_name = 'error';
N=2000
Fr=load('~/simulMCMC/trunk/simuladorBTMM/exp01/conc/conc_ref_0.dat');
fine_m = '~/simulMCMC/trunk/simuladorBTMM/exp01/conc/conc_amostra_';
coarse_m = '~/simulMCMC/trunk/twoStage/simuladorBTMM/exp01/conc/conc_amostra_';

%
Ef=zeros(N,1);
Ec=zeros(N,1);
szp=size(Fr,2)-1;
szt=size(Fr,1)-2;
time=Fr(1:szt,1);
for i=1:N
    namef = [fine_m num2str(i-1) '.dat'];
    namec = [coarse_m num2str(i-1) '.dat'];
    dataf = load(namef);
    datac = load(namec);
    for j=1:szp
        Ff=zeros(szt,1);
        Fc=zeros(szt,1);
        for k=1:szt
            auxf=Fr(k,j)-dataf(k,j);
            auxc=Fr(k,j)-datac(k,j);
            Ff(k)=Ff(k)+auxf*auxf;
            Fc(k)=Fc(k)+auxc*auxc;
        end
%         Ec(i)=Ec(i)+trapz(time,Fc);
%         Ef(i)=Ef(i)+trapz(time,Ff);
        Ec(i)=Ec(i)+sum(Fc);
        Ef(i)=Ef(i)+sum(Ff);
        clear Ff Fc
    end
end
mx=max(max(Ec),max(Ef));
figerro(Ec,Ef,mx);
base=['../figuras/error_' base_name]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r100',base);
%print('-djpeg90',base)
