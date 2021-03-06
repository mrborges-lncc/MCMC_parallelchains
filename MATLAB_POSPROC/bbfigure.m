function bbfigure(N,mu,YMatrix1,nome)
%CREATEFIGURE(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 27-Jul-2018 09:23:08
A = min(min(YMatrix1));
B = max(max(YMatrix1));
%A = min(mu);
if(A>0)
    A=A*.80;
else
    A=1.2*A;
end
if(B>0)
    B=B*1.2;
else
    B=0.8*B;
end
%B = 1.4*max(mu);

% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
    'DataAspectRatio',[1 3*(B-A)/N 1]);
box(axes1,'on');
hold(axes1,'all');

X1 = [1:N]';
% Create multiple lines using matrix input to plot
plot(X1,YMatrix1(:,1,1),'Marker','o','LineStyle','none','Color',[1 0 0]);
plot(X1,mu(1)*ones(1,N),'LineWidth',3,'Color',[0 0 0]);
plot(X1,YMatrix1(:,2,1),'Marker','o','LineStyle','none','Color',[0 1 0]);
plot(X1,mu(2)*ones(1,N),'LineWidth',3,'Color',[0 0 0]);
plot(X1,YMatrix1(:,3,1),'Marker','o','LineStyle','none','Color',[0 0 1]);
plot(X1,mu(3)*ones(1,N),'LineWidth',3,'Color',[0 0 0]);
plot(X1,YMatrix1(:,4,1),'Marker','o','LineStyle','none','Color',[0 0 0]);
plot(X1,mu(4)*ones(1,N),'LineWidth',3,'Color',[0 0 0]);

ylim([A B]);

% Create xlabel
xlabel('Trials','FontSize',18,'FontName','Times New Roman',...
    'FontAngle','italic');

% Create ylabel
ylabel(nome,'FontSize',18,'FontName','Times New Roman',...
    'FontAngle','italic','FontWeight','bold');

