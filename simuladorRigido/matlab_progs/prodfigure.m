function prodfigure(X1, Y1, X2, Y2)
%CREATEFIGURE(X1,Y1,X2,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 16-Mar-2010 12:14:32

% Create figure
figure1 = figure('XVisual',...
    '0x3b (TrueColor, depth 32, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
y_ini=min(Y1);
y_lim=max(Y1)*1.1;
x_ini=min(X1);
x_lim=max(X1);
% y_ini=1.5;
% y_lim=14;

axes1 = axes(...
  'DataAspectRatio',[1/(y_lim-y_ini) 2/(x_lim-x_ini) 1],...
  'Parent',figure1,'FontSize',16,'FontName','Times New Roman');
ylim(axes1,[y_ini y_lim]);
xlim(axes1,[x_ini x_lim]);
box('on');
grid('on');
hold('all');

% Create plot
plot(X1,Y1,'Parent',axes1,'Marker','o','DisplayName','het.');

% Create plot
plot(X2,Y2,'Parent',axes1,'Color',[1 0 0],'DisplayName','hom.');

% Create xlabel
xlabel('$t$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$v_{_{\, Do}}$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.75 0.625 0.07539 0.07296]);

