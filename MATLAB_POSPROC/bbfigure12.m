function bbfigure12(N,Fref1,Fref2,Fref3,xdata,X1,X2,X3,MM,AR,m1,m2,m3)
%CREATEFIGURE(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 27-Jul-2018 09:23:08
% A = min(min(min(X1)));
% A = min(min(min(min(X2))),A);
% A = min(min(min(min(X3))),A);
% B = max(max(max(X1)));
% B = max(max(max(max(X2))),B);
% B = max(max(max(max(X3))),B);
B = max(max(max(Fref2),max(Fref1)),max(Fref3))*1.2;
A = min(min(min(Fref2),min(Fref1)),min(Fref3))-B*0.2;
A = -100;
B = 100;
sz =size(X1);
% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,...
    'FontSize',18,'FontName','Times New Roman',...
    'DataAspectRatio',[1 1*(B-A)/(max(xdata)-min(xdata)) 1]);
box(axes1,'on');
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=[];
f2=[];
f3=[];
for i=MM:sz(1)
    x1=X1(i,:,1);
    x2=X2(i,:,1);
    x3=X3(i,:,1);
    f1=[f1;blackbox(xdata,x1)];
    f2=[f2;blackbox(xdata,x2)];
    f3=[f3;blackbox(xdata,x3)]; 
end
mu1 = mean(f1);
mu2 = mean(f2);
mu3 = mean(f3);
sig1 = sqrt(var(f1));
sig2 = sqrt(var(f2));
sig3 = sqrt(var(f3));
fprintf('\n Desvios médios: sigma_1 = %f; sigma_2 = %f; sigma_3 = %f \n',...
    mean(sig1),mean(sig2),mean(sig3))
% Create multiple lines using matrix input to plot
% for i=0:sz(1)-1
%     x1=X1(N-i,:,1);
%     x2=X2(N-i,:,1);
%     x3=X3(N-i,:,1);
%     Fsn = [blackbox(xdata,x1) blackbox(xdata,x2)...
%         blackbox(xdata,x3)];
%     plot(xdata,Fsn(1:sxd),'ro');
%     plot(xdata,Fsn(sxd+1:2*sxd),'bo');
%     plot(xdata,Fsn(2*sxd+1:3*sxd),'ko');
% end
% Create figure
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-1.1 1.1]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-10 12]);
box(axes1,'on');
hold(axes1,'all');
%
errorbar(xdata,mu1,sig1,'Marker','^','MarkerSize',4,...
    'LineWidth',2,'Color',[1 0 0],'LineStyle','none');
hold on
errorbar(xdata,mu2,sig2,'Marker','s','MarkerSize',4,...
    'LineWidth',2,'Color',[0 0 1],'LineStyle','none');
errorbar(xdata,mu3,sig3,'Marker','o','MarkerSize',4,...
    'LineWidth',2,'Color',[0 0 0],'LineStyle','none');
plot(xdata,Fref1,'r-','LineWidth',3);
plot(xdata,Fref2,'b-','LineWidth',3);
plot(xdata,Fref3,'k-','LineWidth',3);

ylim([A B]);

% Create xlabel
xlabel('x','FontSize',18,'FontName','Times New Roman','FontAngle','italic');

% Create ylabel
ylabel('f(x)','FontSize',18,'FontName','Times New Roman','FontAngle','italic');

% Create textbox
if(AR<99)
    annotation(figure1,'textbox',...
        [0.22 0.15 0.18 0.067],...
        'String',{['AR = ' num2str(AR,'%4.1f') ' %']},...
        'HorizontalAlignment','center',...
        'FontSize',14,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none');
end
if(AR>99)
    sinal = '    ';
    for i=1:length(m1)
        if(m1(i)>=0)
            sinal(i)='+';
        end
    end
    eq = ['$ f(x) = ' num2str(m1(1),'%4.1f') 'x^3'...
        sinal(2) num2str(m1(2),'%4.1f') 'x^2'...
        sinal(3) num2str(m1(3),'%4.1f') 'x'...
        sinal(4) num2str(m1(4),'%4.1f') '$'];
    % Create textbox
    annotation(figure1,'textbox',...
        [0.32 0.82 0.45 0.1],...
        'Interpreter','latex',...
        'String',{eq},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',12,...
        'FontAngle','italic',...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'Color',[1 0 0]);

    % Create textbox
    sinal = '    ';
    for i=1:length(m2)
        if(m2(i)>=0)
            sinal(i)='+';
        end
    end
    eq = ['$f(x) = ' num2str(m2(1),'%4.1f') 'x^3'...
        sinal(2) num2str(m2(2),'%4.1f') 'x^2'...
        sinal(3) num2str(m2(3),'%4.1f') 'x'...
        sinal(4) num2str(m2(4),'%4.1f') '$'];
    annotation(figure1,'textbox',...
        [0.32 0.68 0.45 0.1],...
        'Interpreter','latex',...
        'String',{eq},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',12,...
        'FontAngle','italic',...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'Color',[0 0 1]);

    % Create textbox
    sinal = '    ';
    for i=1:length(m3)
        if(m3(i)>=0)
            sinal(i)='+';
        end
    end
    eq = ['$ f(x) = ' num2str(m3(1),'%4.1f') 'x^3'...
        sinal(2) num2str(m3(2),'%4.1f') 'x^2'...
        sinal(3) num2str(m3(3),'%4.1f') 'x'...
        sinal(4) num2str(m3(4),'%4.1f') '$'];
    annotation(figure1,'textbox',...
        [0.32 0.75 0.45 0.1],...
        'Interpreter','latex',...
        'String',{eq},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontAngle','italic',...
        'FontSize',12,...
        'FontName','Times',...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'Color',[0 0 0]);
end

