% Geracao dos mixings
clear all;
prtprofile = 1;
prtmix     = 2;
prtmix2    = 2;
prtvar     = 1;
prtincl    = 1;
frac = 1;
tic
% [FILENAME, PATHNAME] = uigetfile('/prj/internos/prjmurad/mrborges/rt-penalty-sat-lin/out/*.dat', 'LOAD DATA');
% line_file=sprintf('%s%s', PATHNAME,FILENAME);
% home = '/prj/internos/prjmurad/mrborges/';
home = '/home/mrborges/';
% line_file=[home 'rt-penalty-sat-lin/out/mix03_05_0.dat'];
% line_file=[home 'tracer/out/mean_mix_03_05_0.dat'];
% line_file=[home 'rt-penalty-sat-lin/out/cmix05_0.dat'];
% line_file=[home 'tracer/out/mean_T6_0.dat'];
% line_file=[home 'TRACER_NEW/out/mean_TR23_0.dat'];
line_file=[home 'KL_tracer/out/mix_s1_M10_120x120_0.dat'];
%line_file=[home 'KL_tracer/out/mix_s1_M100_120x120_0.dat'];
%line_file=[home 'KL_tracer/out/mix_s1_M1000_120x120_0.dat'];
%line_file=[home 'KL_tracer/out/mix_s1_M14400_120x120_0.dat'];
line_file='../mix/mix_s05_M40000_400x100_0.dat'
line_file='../mix/mix_e25_a80x20_0.dat'
line_file='../mix/mix_ea80x20_l1_0.dat'
%line_file='../mix/mix_e80x20_l5_0.dat'
%line_file='../mix/mix_e400x100_l5_0.dat'
line_file='../mix/mean_e2000x1000_l50_0.dat';
line_file='../mix/mean_TR1_0.dat'

auxd = 1;
poro = 1*auxd;
dt=10;          % intervalo de tempo
Lx=1000.0;
nx=2000;
dx=Lx/nx;
h=dx*0.5;
NP=10;
inter=1; %inter ==1 => calculo dos mixings interpolando %% cc nao interpola
IA=50;
%
Nrand =150;        % Numero de realizacoes
ini=0;             % Numero da primeira realizacao
fim=ini+Nrand-1;   % Numero da ultima realizacao
satinit = 0.0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut= 75;           % Corte nos tempos finais
cut1=20;           % Corte nos tempos iniciais 05 Lx=2
VY=250;
corte=25;
tracer=11; % ajuste do tempo de simulacao
if(tracer==11)
    corr = (Lx)/(1826.2);
    corr2 = 1;
else
    corr  = 1;
    corr2 = 1;1826.2/Lx;
end
%
for i=size(line_file,2):-1:1
    a=line_file(i);
    if(a=='_')
        break;
    end
end
name_base=line_file(1:i);

time=[];
x=[];
s=[];
v=[];
sat=[];
j=0;
%
for i=ini:fim
    j=j+1;
    file_name = [name_base num2str(i,5) '.dat']
    fn = fopen(file_name,'r');
    test = 1;
    t=0;
    while(test==1)
        textline = fgetl(fn);
        if(textline == -1)
            test=0;
        elseif textline(1:6)=='TEMPO='
            t=t+1;
            if(i==ini)
                time = [time str2num(textline(7:length(textline)))];
            end
            clear dados
            dados = fscanf(fn,'%f');
            dados = [reshape(dados,3,length(dados)/3)'];
            x(:,1)=dados(:,1);
            s(:,t)=dados(:,2);
            v(:,t)=dados(:,3);
        end
    end
    fclose(fn);
    if(i==ini)
        smed=zeros(size(s,1),size(s,2));
        vmed=zeros(size(v,1),size(v,2));
        sat=zeros(Nrand,size(s,1),size(s,2));
    end
    smed = smed+s;
    vmed = vmed+v;
    sat(i+1,:,:) = s;
    clear dados s v
end
smed=smed/Nrand;
vmed=vmed/Nrand;
display('FIM DA LEITURA')
toc

%solucao analitica
vel=(1/poro)*corr*time';
hx=x(2,1,1)-x(1,1,1);
h=0.5*hx;
t=t-cut;
display('TEMPO DE CORTE')
TIME_FIM=time(t)*corr2
TIME_INI=time(cut1)*corr2
clear v
sanl=[];
loc=[];
for i=1:t
    flag=0;
    for j=1:size(x,1)
        if(x(j)+h<=vel(i))
            sanl(j,i)=1.0;
        else
            sanl(j,i)=0.0;
            if(flag==0)
                loc(i)=j;
                flag=1;
            end
        end
    end
end

%%% calculo da integral para determinar o mixing
fat=2*sqrt(pi);
mix=zeros(1,t);
if(inter==1)
%%% calculo da integral para determinar o mixing com interpolacao
    xx=[min(x):hx/(IA):max(x)];
    %
    for i=1:t
        yy=interp1(x,smed(:,i),xx,'spline');
        for j=1:size(xx,2)
            if(vel(i)<=xx(j))
                loc = j;
                break
            end
        end
        mix(i) = trapz(xx(loc:end),yy(loc:end));
    end
    if(mix(1,max(size(mix))))==0
        i
        break
    end
  mixi=zeros(Nrand,t);
  for k=1:Nrand
    for i=1:t
        yy=interp1(x,sat(k,:,i),xx,'spline');
        for j=1:size(xx,2)
            if(vel(i)<=xx(j))
                loc = j;
                break
            end
        end
        mixi(k,i) = trapz(xx(loc:end),yy(loc:end))';
    end
  end
else
%%% calculo da integral para determinar o mixing sem interpolacao
    for i=1:t
        mix(i) = trapz(x(loc(i):end),smed(loc(i):end,i));
    end
end

mix=fat*mix*auxd;
mixi=fat*mixi*auxd;
clear sat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('FIM DO CALCULO DOS MIXINGS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Estatisticas')
%vY = log(mixi(:,int16(size(mixi,2)*.75)));
vY = log(mixi(1,2:end));
mvY = mean(vY);
vvY = var(vY);
%hist(vY);
smix = sqrt((vvY));
xmix = mvY;
ic1 = exp(xmix+1.*smix);
ic2 = exp(xmix-1.*smix);
% ic1 = exp(xmix+1.96*smix/sqrt(Nrand));
% ic2 = exp(xmix-1.96*smix/sqrt(Nrand));
% break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(prtprofile==1)
    % Create figure PROFILES
    figure1 = figure();

    % Create axes
    varx=max(x+h+hx)*1.1;
    vary=(max(max(smed))*1.3);
    axes1 = axes('Parent',figure1,...
        'FontSize',16,...
        'FontName','Times',...
        'DataAspectRatio',[varx*0.5 vary 0.01]);
    ylim(axes1,[0.0 vary]);
    xlim(axes1,[0.0 varx]);

    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim([0 1.2]);
    box('on');
    hold('all');
    for T=1:corte:t
    plot(x,smed(:,T),'o','MarkerSize',4,'Parent',axes1)
    hold on
    plot(x,sanl(:,T),'k','LineStyle','--','Parent',axes1)
    end
    % Create xlabel
    xlabel('$x_1$','Interpreter','latex','FontSize',18,'FontName','Times');

    % Create ylabel
    ylabel('$\langle S \rangle$','Interpreter','latex','FontSize',18,...
        'FontName','Times');
    %% Create legend
    legend1 = legend(...
      axes1,{...
      'computed','hom.'},...
      'FontName','times',...
      'FontSize',14,...
      'Position',[0.3 0.635 1.0 0.1849],...
      'Interpreter','latex');
    legend(legend1,'boxoff');
    sinal=0;
    for i=size(name_base,2):-1:1
        str=name_base(i);
        if(str=='/')
            sinal=sinal+1;
        end
        if(sinal==2)
            break;
        end
    end
    name_fig = [name_base(1:i) 'figuras/mixing/prof_' name_base(i+5:end) num2str(Nrand,5)]

    % print('-djpeg99',name_fig);
    print('-depsc',name_fig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%% Create figure VARIANCIAS
if(prtvar==1)
    figure2 = figure();

    % Create axes
    varx=max(x+h+hx)*1.1;
    vary=(max(max(vmed))*1.4);
    axes2 = axes('Parent',figure2,...
        'FontSize',16,...
        'FontName','Times',...
        'DataAspectRatio',[varx*0.5 vary 0.01]);
    ylim(axes2,[0.0 vary]);
    xlim(axes2,[0.0 varx]);

    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim([0 1.2]);
    box('on');
    hold('all');
    for T=1:corte:t
    plot(x,vmed(:,T),'o','MarkerSize',4,'Parent',axes2)
    hold on
    c=max(max(vmed));
    plot(x,c*1.1*sanl(:,T),'k','LineStyle','--','Parent',axes2)
    end
    % Create xlabel
    xlabel('$x_1$','Interpreter','latex','FontSize',18,'FontName','Times');

    % Create ylabel
    ylabel('$\sigma^{2}_{S}$','Interpreter','latex','FontSize',18,...
        'FontName','Times');
    %% Create legend
    legend2 = legend(...
      axes2,{...
      'computed','hom.'},...
      'FontName','times',...
      'FontSize',14,...
      'Position',[0.3 0.635 1.0 0.1849],...
      'Interpreter','latex');
    legend(legend2,'boxoff');

    name_fig = [name_base(1:i) 'figuras/mixing/var_' name_base(i+5:end) num2str(Nrand,5)]

    % print('-djpeg99',name_fig);
    print('-depsc',name_fig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIGURAS MIXINGS
    clear x y
    time=corr2*time;
    x=mix';
    y=time(1:end-cut)';
    clear x y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=(time(cut1+1:1:end-cut)');
    X=[];
    Y=[];
    k=0;
    % for i=0:Nrand-1
    %     y=(mixi(i+1,cut1+1:1:end)');
    %     for j=1:size(y,1)
    %         y=(mixi(i+1,cut1+1:1:end)');
    %         k=k+1;
    %         X(k)=x(j);
    %         Y(k)=y(j);
    %     end
    % end
    y=(mixi(1,cut1+1:1:end)');
    for j=1:size(y,1)
        for i=1:Nrand
            y=(mixi(i,cut1+j)');
            k=k+1;
            X(k)=x(j);
            Y(k)=y;
        end
    end
    clear x y
    aux=corrcoef([X' Y']);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x=log(time(cut1+1:1:end-cut)');
    % y=log(mix(cut1+1:1:end)');
    % %[P,S,MU]=polyfit(x,y,1);
    % P=polyfit(x,y,1);
    % aux=corrcoef([x y]);
    % clear x y
    % x = 0:max(time)/500:max(time);
    % y = exp(P(2))*x.^(P(1));
    % equacao = ['$\ell(t) = ' num2str(exp(P(2)),2) 't^{' num2str(P(1),2) '}'];
    % equacao=[equacao '\quad R^{2} = ' num2str(aux(1,2)*aux(1,2),3) '$'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=(time(cut1+1:1:end-cut)');
    y=(mix(cut1+1:1:end)');
    %[P,S,MU]=polyfit(x,y,1);
    P=polyfit(x,y,2);
    aux=corrcoef([x y]);
    clear x y
    x = 0:max(time)/500:max(time);
    y = (P(1))*x.^2+(P(2))*x+P(3);
    equacao = ['$\ell(t) = ' num2str((P(1)),2) 't^{2}+' num2str(P(2),2) 't' num2str(P(3),2)];
    equacao=[equacao '\quad R^{2} = ' num2str(aux(1,2)*aux(1,2),3) '$'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create figure
if(prtmix==2)
    if(frac==1)
        mix_frac(time(cut1+1:end-cut),mix(cut1+1:end),time,cut,VY);
    else
        mix_regr(time(cut1+1:end-cut),mix(cut1+1:end),time,cut,VY);
    end
    %%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sinal=0;
    for i=size(name_base,2):-1:1
        str=name_base(i);
        if(str=='/')
            sinal=sinal+1;
        end
        if(sinal==2)
            break;
        end
    end

    name_fig = [name_base(1:i) 'figuras/mixing/lt_' name_base(i+5:end) num2str(Nrand,5)]

    % print('-djpeg99',name_fig);
    print('-depsc',name_fig);
end
if(prtmix==1)
    figure3 = figure();

    % Create axes
    varx=max(time(1:end-cut))*1.1;
    %vary=max(mix)*1.3;
    vary=VY;
    axes3 = axes('Parent',figure3,...
        'FontSize',16,...
        'FontName','Times',...
        'DataAspectRatio',[varx*0.5 vary 0.01]);
    ylim(axes3,[0.0 vary]);
    xlim(axes3,[0.0 varx]);

    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim([0 1.2]);
    box('on');
    hold('all');
    plot(time(1:end-cut),mix','o','MarkerSize',4,'Color',[0 0 0],'Parent',axes3);
    %plot(time(1:end-cut),mix','-','LineWidth',2,'Color',[0 0 0],'Parent',axes3);
    plot(x,y,'r','LineStyle','-','LineWidth',2,'Parent',axes3);
    plot(time(1:end-cut),ic1,'LineWidth',1,'LineStyle','--','Color',[0.8471 0.1608 0],'Parent',axes3);
    plot(time(1:end-cut),ic2,'LineWidth',1,'LineStyle','-.','Color',[0.8471 0.1608 0],'Parent',axes3);
    x5=x;
    y5=y;
    % Create xlabel
    xlabel('$t$','Interpreter','latex','FontSize',18,'FontName','Times');

    % Create ylabel
    ylabel('$\ell\, (t)$','FontSize',18,...
        'FontName','Times','Interpreter','latex');
    %% Create legend
    legend3 = legend(...
      axes3,'String',{...
      '$\mbox{computed}$',equacao,'$$\mu+\sigma$$','$$\mu-\sigma$$'},...
      'FontName','times',...
      'FontSize',14,...
      'Position',[-0.05 0.565 1.0 0.2],'Interpreter','latex');
    %legend(legend3,'interpreter','latex');
    legend(legend3,'boxoff');
    %%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sinal=0;
    for i=size(name_base,2):-1:1
        str=name_base(i);
        if(str=='/')
            sinal=sinal+1;
        end
        if(sinal==2)
            break;
        end
    end

    name_fig = [name_base(1:i) 'figuras/mixing/lt_' name_base(i+5:end) num2str(Nrand,5)]

    % print('-djpeg99',name_fig);
    print('-depsc',name_fig);
end
%%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sinal=0;
for i=size(name_base,2):-1:1
    str=name_base(i);
    if(str=='/')
        sinal=sinal+1;
    end
    if(sinal==2)
        break;
    end
end
%%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome_out=[name_base(1:i) 'out/lt_' name_base(i+5:end) num2str(Nrand,5) '.dat']
fout = fopen(nome_out,'wt');

for j=1:size(X,2)%t
%    fprintf(fout,'%g %g\n',time(j),mix(j));
    fprintf(fout,'%g %g\n',X(j),Y(j));
end
fclose(fout);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% verificar a inclinacao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(prtincl==1)
    gamma=0.75
    NV=int16(size(mix,2)/NP-1e-7);
    incl=[];
    xinc=[];
    for k=1:NV-1
        clear x y
        idi=NP*(k-1)+1;
        idf=NP*k;
        if(k==1)
            idi=idi+1;
        end
        x=log(time(idi:idf)');
        y=log(mix(idi:idf)');
        P=polyfit(x,y,1);
        incl(k)=P(1);
        xinc(k)=time(int16(0.5*NP+idi));
    end
    figure10 = figure();
    varx=max(xinc)*1.1;
    vary=max(incl)*1.1;
    varyy=min(incl)*.9;
    varxx=min(xinc)*.9;
    vary=1.1;
    varyy=.4;

    axes10 = axes('Parent',figure10,...
        'FontSize',16,...
        'FontName','Times',...
        'DataAspectRatio',[(varx-varxx)*0.5 (vary-varyy) 0.01]);
    ylim(axes10,[varyy vary]);
    xlim(axes10,[varxx varx]);
    z=xinc*0.0+gamma;
    box('on');
    hold('all');
    plot(xinc,z,'-','MarkerSize',4,'Parent',axes10);
    plot(xinc,incl,'o','MarkerSize',6,'Parent',axes10);
    xlabel('$t_c$','Interpreter','latex','FontSize',20,'FontName','Times');
    grid('on');

    % Create ylabel
    ylabel('$\gamma$','Interpreter','latex','FontSize',20,...
        'FontName','Times');
    name_fig = [name_base(1:i) 'figuras/mixing/incl_' name_base(i+5:end) num2str(Nrand,5)]

    % print('-djpeg99',name_fig);
    print('-depsc',name_fig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURAS MIXINGS II
if(prtmix2==1)

    jump=4
    x=[];
    y=[];
    for i=1:max(size(X))/jump
        k = (i)*jump;
        x(i) = X(k);
        y(i) = Y(k);
    end
    for k=1:1%Nrand

        if(k==1)
            time=corr2*time;
            % Create figure
           figure5 = figure();

            % Create axes
            varx=max(time(1:end-cut))*1.1;
            %vary=max(mix)*1.3;
            vary=VY;
            axes5 = axes('Parent',figure5,...
                'FontSize',16,...
                'FontName','Times',...
                'DataAspectRatio',[varx*0.5 vary 0.01]);
            ylim(axes5,[0.0 vary]);
            xlim(axes5,[0.0 varx]);

            % Uncomment the following line to preserve the Y-limits of the axes
            % ylim([0 1.2]);
            box('on');
            hold('all');
        end
        plot(x,y,'o','MarkerSize',4,'Color',[0 0 0],'Parent',axes5);
        hold on
        plot(time(1:end-cut),mix','-','LineWidth',2,'Color',[1 0 0],'MarkerSize',4,'Parent',axes5);
        plot(time(1:end-cut),ic1,'LineWidth',1,'LineStyle','--','Color',[1 1 0],'Parent',axes5);
        plot(time(1:end-cut),ic2,'LineWidth',1,'LineStyle','-.','Color',[1 1 0],'Parent',axes5);
        % plot(time(1:end-cut),mix'+desv','LineWidth',1,'LineStyle','--','Color',[0.8471 0.1608 0],'Parent',axes5);
        % plot(time(1:end-cut),mix'-desv','LineWidth',1,'LineStyle','--','Color',[0.8471 0.1608 0],'Parent',axes5);
        % Create xlabel
    %    if(k==Nrand)
        if(k==1)
            %plot(x5,y5,'k','LineStyle','-','LineWidth',2,'Parent',axes5);
            %xlabel('$t$','Interpreter','latex','FontSize',18,'FontName','Times');
            % Create xlabel
            xlabel('$t$','Interpreter','latex','FontSize',18,'FontName','Times');

            % Create ylabel
            ylabel('$\ell(t)$','Interpreter','latex','FontSize',18,...
                'FontName','Times');
            %% Create legend
            legend5 = legend(...
              axes5,{...
              '$computed$',equacao,'$$\mu+\sigma$$','$$\mu-\sigma$$'},...
              'FontName','times',...
              'FontSize',14,...
              'Position',[0.05 0.6 1.0 0.1849],...
              'Interpreter','latex');
            legend(legend5,'boxoff');

    %%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sinal=0;
    for i=size(name_base,2):-1:1
        str=name_base(i);
        if(str=='/')
            sinal=sinal+1;
        end
        if(sinal==2)
            break;
        end
    end
        name_fig = [name_base(1:i) 'figuras/mixing/lt_all_' name_base(i+5:end) num2str(Nrand,5)]

        % print('-djpeg99',name_fig);
        print('-depsc',name_fig);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

if(prtmix2==2)
    jump=4
    x=[];
    y=[];
    for i=1:max(size(X))/jump
        k = (i)*jump;
        x(i) = X(k);
        y(i) = Y(k);
    end
    if(frac==1)
        mix_frac(X,Y,time,cut,VY)
    else
        mix_regr(X,Y,time,cut,VY)
    end
    sinal=0;
    for i=size(name_base,2):-1:1
        str=name_base(i);
        if(str=='/')
            sinal=sinal+1;
        end
        if(sinal==2)
            break;
        end
    end
        name_fig = [name_base(1:i) 'figuras/mixing/lt_all_' name_base(i+5:end) num2str(Nrand,5)]

        % print('-djpeg99',name_fig);
        print('-depsc',name_fig);
end
clear;