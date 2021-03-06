function twoEphiplot(X1, X2, Y1, Y2, xname, yname, dat1, dat2, E0, c, n)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data
n1 = num2str(E0,'%1.2e');
n1 = [n1(1:length(n1)-4) '\times 10^{' n1(end-2:end) '}'];
eq = ['$' yname '(' xname ') = ' n1  '\displaystyle \exp\left(' num2str(c,'%1.2f') xname '\right)$'];
A = min(X1); B = max(X1);
C = min(min(Y1)); D = max(max(Y1));
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
plot1 = plot(X1,Y1,'DisplayName',dat1,'MarkerSize',8,'Marker','o',...
    'LineWidth',2,...
    'LineStyle','none',...
    'Color',[0 0.45 0.74]);
plot(X2,Y2,'DisplayName',eq,'LineWidth',3,...
    'Color',[0.93 0.69 0.125]);

% Create ylabel
yname = ['$' yname '$'];
ylabel({yname},'FontSize',22,'Interpreter','latex');

% Create xlabel
xname = ['$' xname '$'];
xlabel({xname},'FontSize',22,'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',asp,...
    'FontName','Times','FontSize',14,'XMinorTick','on','YMinorTick',...
    'on','ZMinorTick','on');
xlim([A B]);
ylim([C D]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',12,'FontName','Times','EdgeColor','none',...
    'Color','none','Interpreter','latex','Location','NorthEast');


