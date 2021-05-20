function sat_figure_corte(xdata1,zdata1,cc,dd,A,B)
%CREATEFIGURE(XDATA1,YDATA1,ZDATA1)
%  XDATA1:  surface xdata
%  YDATA1:  surface ydata
%  ZDATA1:  surface zdata

%  Auto-generated by MATLAB on 16-Feb-2011 12:02:36

% Create figure
figure1 = figure('PaperType','<custom>','PaperSize',[18.92 14.57],...
    'PaperOrientation','landscape',...
    'PaperUnits','centimeters');
if(abs(A-B)<1e-5)
    posi=[0.1 0.11 0.775 0.815];
else
    posi=[0.1 0.11 0.775 0.815];
end
% Create axes
axes1 = axes('Parent',figure1,'YTick',[0 dd*0.2 dd*0.4 dd*0.6 dd*0.8 dd*1],...
    'YMinorTick','on',...
    'YGrid','on',...
    'XTick',[0 A*0.2 A*0.4 A*0.6 A*0.8 A*1],...
    'XMinorTick','on',...
    'XGrid','on',...
    'FontSize',18,...
    'FontName','times',...
    'DataAspectRatio',[1 2*abs(dd-cc)/A 1],...
    'CLim',[cc dd],...
    'Position',posi);
% Create axes
xlim(axes1,[0 A]);
dd=dd*1.1;
ylim(axes1,[cc dd]);
grid(axes1,'on');
hold(axes1,'all');
box on

% Create surf
plot(xdata1,zdata1,...
    'MarkerEdgeColor',[0 0 0],'Marker','o','LineWidth',1,'Color',[0 0 1],'LineStyle','none','DisplayName','Numerical');
    hold on
x=[0:0.0001:1.0];    
y=0.0*x;
for i=1:length(y)
    if(x(i)>=0.85)
        if(x(i)<=0.95)
            y(i)=1.0;
        end
    end
end
plot(x,y,...
    'MarkerEdgeColor',[0 0 0],'Marker','none','LineWidth',1,'Color',[0 0 1],'DisplayName','Exact');
% Create zlabel
ylabel('s','Visible','on','FontSize',24,'FontName','times','Interpreter','tex');
xlabel('\textit{x}','Visible','on','FontSize',24,'FontName','times','Interpreter','latex');
view(0,90.0)
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.155 0.64 0.1635 0.1412]);

