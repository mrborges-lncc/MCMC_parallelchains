function sat_figure_vista(xdata1,ydata1,zdata1,A,B,aa,bb,cc,dd,C,D,t)
%CREATEFIGURE(XDATA1,YDATA1,ZDATA1)
%  XDATA1:  surface xdata
%  YDATA1:  surface ydata
%  ZDATA1:  surface zdata

%  Auto-generated by MATLAB on 16-Feb-2011 12:02:36

% Create figure
figure1 = figure('PaperType','<custom>','PaperSize',[18.92 14.57],...
    'PaperOrientation','landscape',...
    'PaperUnits','centimeters');

% Create axes
axes1 = axes('Parent',figure1,'ZTick',[cc 0.5*(dd+cc)  dd],...
    'YTick',[0 B*0.2 B*0.4 B*0.6 B*0.8 B*1],...
    'YMinorTick','on',...
    'YGrid','on',...
    'XTick',[0 A*0.2 A*0.4 A*0.6 A*0.8 A*1],...
    'XMinorTick','on',...
    'XGrid','on',...
    'FontSize',16,...
    'FontName','times',...
    'DataAspectRatio',[1 1 4*abs(D-C)/A],...
    'CLim',[0.2 0.85]);
% Create axes
xlim(axes1,[0 A]);
ylim(axes1,[0 B]);
grid(axes1,'on');
hold(axes1,'all');
box on

% Create surf
contour3(xdata1,ydata1,zdata1,8);%,...
    %'Parent',axes1,'LineWidth',1.0);%,'ShowText','on');
hold on
surf(xdata1,ydata1,zdata1);
    shading('interp');
    alpha(0.5);
    hold off
    view(55,23);
% Create zlabel
zlabel('S_{w}','Visible','off','FontSize',20,'FontName','times');

% Create colorbar
%colorbar('peer',axes1,[0.8522 0.1071 0.04216 0.8071],'FontSize',14);
%
% Create textbox
tempo = num2str(t,'%1.4f');
if(abs(A-B)<1e-8)
    annotation(figure1,'textbox',[0.4331 0.9238 0.1151 0.03881],...
        'String',{['t =' tempo]},...
        'HorizontalAlignment','center',...
        'FitBoxToText','off',...
        'LineStyle','none');
else
    annotation(figure1,'textbox',[0.45 0.775 0.1151 0.03881],...
        'String',{['t =' tempo]},...
        'HorizontalAlignment','center',...
        'FitBoxToText','off',...
        'LineStyle','none');
end



