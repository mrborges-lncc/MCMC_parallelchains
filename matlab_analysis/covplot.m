function covplot(X1, X2, Y1, Y2, xname, yname, dat1,dat2, n)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data
A = min(min(X1),min(X2)); B = max(max(X1),max(X2));
C = min(min(Y1),min(Y2)); D = round(max(max(Y1),max(Y2)))*1.1;
if abs(D)<1e-18
    D = 0.001;
else
    if abs(D-C)<1e-18
        C = C - D*0.1;
        D = D + D*0.1;
    end
end
asp = [(B-A)/(D-C) 1 1];%  Auto-generated by MATLAB on 26-Feb-2021 10:45:07

% Create figure
figure1 = figure(n);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X2,Y2,'DisplayName',dat2,'LineWidth',3,...
    'Color',[0.93 0.69 0.125]);
plot(X1,Y1,'DisplayName',dat1,'MarkerSize',6,'Marker','o',...
    'LineWidth',1,...
    'LineStyle','none',...
    'Color',[0 0.45 0.74]);
% Create ylabel
ylabel({yname},'FontSize',24,'Interpreter','latex');

% Create xlabel
xlabel({xname},'FontSize',24,'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',asp,'TickLabelInterpreter','latex',...
    'FontName','Times','FontSize',16,'XMinorTick','on','YMinorTick',...
    'on','ZMinorTick','on');
xlim([A B]);
ylim([C D]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',12,'FontName','Times','EdgeColor','none',...
    'Color','none','Interpreter','latex','Location','NorthEast');

