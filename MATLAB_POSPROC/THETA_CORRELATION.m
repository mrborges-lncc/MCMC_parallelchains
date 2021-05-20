%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELACAO ENTRE THETAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FILENAME, PATHNAME] =uigetfile('../twoStage/select_thetas/*.dat',...
    'LOAD DATA');
nc=getNAME(FILENAME)
namet1 = FILENAME(1:nc);
[FILENAME, PATHNAME] =uigetfile('../twoStage/select_thetas/*.dat',...
    'LOAD DATA');
namet2 = FILENAME(1:nc);
ini=00;
fim=21;
Nthetas1=50;
Nthetas2=50;
Nthetas=min(Nthetas1,Nthetas2);
NT =fim-ini+1;
cut1=0;
cut2=40;
NTth=Nthetas-cut1-cut2;
base_name='KL20'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ELEMENTOS COMPARADOS')
1+cut1:Nthetas-cut2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat1=zeros(1,Nthetas1);
dat2=zeros(1,Nthetas2);
data1=zeros(NT,NTth);
data2=zeros(NT,NTth);
for i=1:NT
    name = [PATHNAME namet1 num2str(ini+i-1) '.dat'];
    dat1(1,:)=load(name);
    data1(i,:)=dat1(1+cut1:Nthetas-cut2);
    name = [PATHNAME namet2 num2str(ini+i-1) '.dat'];
    dat2(1,:)=load(name);
    data2(i,:)=dat2(1+cut1:Nthetas-cut2);
end
mu1 = (mean(data1,2));
mu2 = (mean(data2,2));
var1= var(data1');
var2= var(data2');
x=[1:1:NT];
figMEANS(x,[mu1 mu2],0.0);
base=['../figuras/Meanthetas_' base_name];
% set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
figVARS(x,[var1' var2'],1.0);
base=['../figuras/Varthetas_' base_name];
% set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cov=zeros(NT,1);
for i=1:NT
    sum1=0.0;
    mu1=mean(data1(i,:));
    mu2=mean(data2(i,:));
    for j=1:NTth
        sum1 = sum1+(data1(i,j)-mu1)*(data2(i,j)-mu2);
    end
    cov(i)=sum1/NTth;
end
MAX1=max(cov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cov=zeros(NTth,1);
for i=1:NTth
    sum1=0.0;
    mu1=mean(data1(:,i));
    mu2=mean(data2(:,i));
    for j=1:NT
        sum1 = sum1+(data1(j,i)-mu1)*(data2(j,i)-mu2);
    end
    cov(i)=sum1/NT;
end
MAX2=max(cov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vet=zeros(NT*NTth,2);
k=1;
for i=1:NT
    for j=1:NTth
        vet(k,1)=(data1(i,j));
        vet(k,2)=(data2(i,j));
        k=k+1;
    end
end
cr=corrcoef(vet(:,1),vet(:,2));
figTHETAS(vet(:,1),vet(:,2),cr);
base=['../figuras/thetas_' base_name];
% set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);

    