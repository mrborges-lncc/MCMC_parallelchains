clear;
base_name = 'prodF';
file_name = ['../MCMC/error/niveis_' base_name '.dat']
dados=load(file_name);
sz=size(dados);

selN=[];
noselN=[];
for i=1:sz(1)
    if(dados(i,2)==1)
        selN=[selN; dados(i,1)];
    else
        noselN=[noselN; dados(i,1)];
    end
end
disp('Tamanho de selN')
NA=size(selN,1)
N=max(dados(:,1))
ma=NA;
mi=0;
% 
freq=zeros(N+1,1);
% 
for i=0:N
    for j=1:NA
        if(selN(j)==i)
            freq(i+1)=freq(i+1)+1;
        end
    end
end
figure1=figure(1);
%axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman')%,...
    %'DataAspectRatio',[max(freq(:,1))*2 N 2]);
% Create axes
x=[-1:1:N];
axes1 = axes('Parent',figure1,...
    'XTick',x,...
    'FontSize',16,...
    'FontName','Times New Roman',...
    'CLim',[1 2]);
box(axes1,'on');
hold(axes1,'all');
xlim(axes1,[-0.5 N+0.5])
ylim(axes1,[0 max(freq(:,1))])
% 
%x=[0:1:N];
%plot(x,freq(:,1));
%
edges=[0:0.5:20];
hist(selN,edges)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','none','EdgeColor','k')% 
% Create xlabel
xlabel('Level','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman','FontAngle','italic');

% Create ylabel
ylabel('Freq.','FontWeight','bold','FontSize',20,...
    'FontName','Times New Roman','FontAngle','italic');

% Print
base=['../figuras/niveis_' base_name];
%print('-djpeg90',base)
print('-depsc','-r100',base)
