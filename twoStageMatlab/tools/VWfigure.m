function VWfigure(X1, Y1, Y2, Y3, V1, V2, V3)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 24-Oct-2018 15:46:25
B   = max(max(max(V1),max(V2)),max(V3));
B   = max(max(max(max(Y1),max(Y2)),max(Y3)),B);
aux = B * 0.1;
B   = B + aux;
A   = min(min(min(V1),min(V2)),min(V3));
A   = min(min(min(min(Y1),min(Y2)),min(Y3)),A) - aux;
aux = max(abs(B-1),abs(A-1));
% A = 0.8;%1-aux;
% B = 1.3;%1+aux;
C = 0.0;%min(X1);
D = max(X1);
% Create figure
figure1 = figure('OuterPosition',[95 306 481 351]);

% Create axes
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'FontSize',16,'FontName','Times',...
    'DataAspectRatio',[((D-C)/(B-A)) 2 1]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
for i=1:size(Y1,2)
    plot(X1,Y1(:,i),'LineWidth',3,'Color',[0 0 0],'LineStyle',':',...
        'MarkerSize',6,'Marker','none','DisplayName','$\mathsf{W}_1^{1/2}$');
end
for i=1:size(V1,2)
    plot(X1,V1(:,i),'LineWidth',3,'Color',[0 0 0],'LineStyle',...
        '--','DisplayName','$\mathsf{V}_1^{1/2}$');
end
% for i=1:size(Y2,2)
%     plot(X1,Y2(:,i),'LineWidth',3,'Color',[1 0 0],'LineStyle','-','DisplayName','$C_2$');
% end
% for i=1:size(Y3,2)
%     plot(X1,Y3(:,i),'LineWidth',3,'Color',[0 0 1],'LineStyle','-','DisplayName','$C_3$');
% end
% plot(X1,1+0.0*Y1(:,i),'LineWidth',1,'LineStyle',':','Color',[0.5 0.5 0.5],'DisplayName','$1.0$');
% plot(X1,1.2+0.0*Y1(:,i),'LineWidth',1,'LineStyle','--','Color',[0.5 0.5 0.5],'DisplayName','$1.2$');

% Create xlabel
xlabel('Interaction number','Interpreter','none','FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times');

% Create ylabel
ylabel('Estimative','Interpreter','latex','FontSize',20,'FontName','Times');

ylim([A B])
xlim([C D])

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex','FontSize',14,'EdgeColor','none',...
    'Color','none');
