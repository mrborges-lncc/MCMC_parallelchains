function sat_figure(xdata1,ydata1,zdata1,A,B,aa,bb,cc,dd,C,D,t,clim1,clim2)
%CREATEFIGURE(XDATA1,YDATA1,ZDATA1)
%  XDATA1:  surface xdata
%  YDATA1:  surface ydata
%  ZDATA1:  surface zdata

%  Auto-generated by MATLAB on 16-Feb-2011 12:02:36

% Create figure
figure1 = figure('PaperType','<custom>','PaperSize',[18.92 10.57],...
    'PaperOrientation','landscape',...
    'PaperUnits','centimeters');
if(abs(A-B)<1e-5)
    posi=[0.1 0.11 0.775 0.815];
else
    if(abs(A/B)<2.1)
        posi=[0.1 0.23 0.775 0.815];
    else
        posi=[0.1 0.35 0.775 0.815];
    end
end

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
    'CLim',[clim1 clim2],...
    'Position',posi);
% Create axes
xlim(axes1,[0 A]);
ylim(axes1,[0 B]);
grid(axes1,'on');
%hold(axes1,'all');
box on

% Create surf
    surf(xdata1,ydata1,zdata1);
    shading('interp');
%    shading('flat');
%    shading('faceted');
%    alpha(0.5); 
    hold on
%    contour(xdata1,ydata1,zdata1,8,...
%    'Parent',axes1,'LineWidth',1.0);%,'ShowText','on'); 
    hold off
% Create zlabel
zlabel('S_w','Visible','off','FontSize',20,'FontName','times');
view(0,90.0)
% Create colorbar
%if(abs(A-B)<1e-5)
%    colorbar('peer',axes1,[0.8522 0.1071 0.04216 0.8071],'FontSize',15,'FontName','times');
%else
%    if(abs(A/B)<2.1)
%        colorbar('peer',axes1,[0.9 0.377 0.04216 0.52],'FontSize',15,'FontName','times');
%    else
%        colorbar('peer',axes1,[0.9 0.632 0.04216 0.26],'FontSize',15,'FontName','times');
%    end
%end
%
% % Create textbox
% tempo = num2str(t,'%1.3f');
% if(abs(A-B)<1e-5)
%     annotation(figure1,'textbox',[0.4331 0.9238 0.1151 0.03881],...
%         'String',{['t =' tempo]},...
%         'HorizontalAlignment','center',...
%         'FitBoxToText','off',...
%         'LineStyle','none','FontName','times','FontSize',12);
% else
%     if(abs(A-B)<2.1)
%         annotation(figure1,'textbox',[0.45 0.92 0.1151 0.03881],...
%             'String',{['t =' tempo]},...
%             'HorizontalAlignment','center',...
%             'FitBoxToText','off',...
%             'LineStyle','none','FontName','times','FontSize',12);
%     else
%         annotation(figure1,'textbox',[0.45 0.92 0.1151 0.03881],...
%             'String',{['t =' tempo]},...
%             'HorizontalAlignment','center',...
%             'FitBoxToText','off',...
%             'LineStyle','none','FontName','times','FontSize',12);
%     end        
% end



