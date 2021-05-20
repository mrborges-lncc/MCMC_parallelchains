function fig_conc(X1,YMatrix1,str)
%CREATEFIGURE(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 09-Nov-2011 13:47:52

% Create figure
figure1 = figure('XVisual',...
    '0x23 (TrueColor, depth 32, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
axes1 = axes('Parent',figure1,'FontSize',12,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Marker','o','Color',[0 0 0],...
    'Parent',axes1,'LineStyle','none',...
    'DisplayName',str);

% Create xlabel
xlabel('time','Interpreter','latex','FontSize',14,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('Conc.','Interpreter','latex','FontSize',14,...
    'FontName','Times New Roman');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[0 0 0],...
    'Position',[0.7382 0.748 0.1073 0.128],...
    'FontAngle','italic','Color','none');

