function two2Dplot(X1, Y1, xname, yname, dat1,dat2, n)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data
A = min(X1); B = max(X1);
C = min(min(Y1)); D = max(max(Y1));
if abs(D)<1e-6
    D = 0.001;
else
    if abs(D-C)<1e-6
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
plot1 = plot(X1,Y1,'LineWidth',2);
set(plot1(1),'DisplayName',dat1,...
    'Color',[0 0.45 0.74]);
set(plot1(2),'DisplayName',dat2,...
    'Color',[0.85 0.33 0.1]);

% Create ylabel
ylabel({yname},'FontSize',18,'Interpreter','latex');

% Create xlabel
xlabel({xname},'FontSize',18,'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',asp,...
    'FontName','Times','FontSize',12,'XMinorTick','on','YMinorTick',...
    'on','ZMinorTick','on');
xlim([A B]);
ylim([C D]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',12,'FontName','Times','EdgeColor','none',...
    'Color','none','Interpreter','latex');

