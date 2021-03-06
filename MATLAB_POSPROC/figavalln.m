function figavalln(M, Y1, Y2, st)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 19-Dec-2011 08:09:49
X1=[1:M]';
Y1=Y1(1:M);
Y2=Y2(1:M);
B=max(Y1)*1.1;
A=0.0;%min(Y1);
%B=1.5e+1;
C=min(X1);
D=max(X1);

% Create figure
figure1 = figure('XVisual',...
    '0x27 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
dasp=[1 (B-A)/(D-C) 200];
%axes1 = axes('Parent',figure1,'YScale','linear','YMinorTick','on',...
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');
xlim(axes1,[0 D]);
ylim(axes1,[A B]);

% Create semilogy
% plot(X1,Y1,'Parent',axes1,'MarkerSize',8,'Marker','o','LineWidth',2,...
%     'DisplayName','Independent sampling',...
%     'Color',[0 0 0]);
semilogy(X1,Y1,'Parent',axes1,'MarkerSize',4,'Marker','none','LineWidth',2,...
    'DisplayName','Exponential',...
    'Color',[1 0 0]);
hold on
semilogy(X1,Y2,'Parent',axes1,'MarkerSize',4,'Marker','none','LineWidth',2,...
    'DisplayName','Square Exp.',...
    'Color',[0 0 1]);
% Create xlabel
xlabel('$M$','Interpreter','latex','FontSize',20,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$\ln( \lambda\ )$','Interpreter','latex','FontSize',20,...
    'FontName','Times New Roman');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.631 0.835 0.07064 0.05092],'Box','off');

% Create textbox
% nome = ['\Phi_{' num2str(M,5) '} = ' num2str(st,'%8.2f')];
nome = ['G({' num2str(M,5) '}) = ' num2str(st,'%8.2f')];
% annotation(figure1,'textbox',[0.6135 0.6929 0.1972 0.0881],...
%     'String',{nome},...
%     'FontSize',12,...
%     'FontName','Times New Roman',...
%     'FontAngle','italic',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

