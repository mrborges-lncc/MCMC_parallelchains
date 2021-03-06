function figPeaceman(x,y,xx,yy)
    xmin = min(min(x),min(xx));
    xmax = max(max(xx),max(x));
    ymin = min(min(yy),min(y));
    ymax = max(max(yy),max(y));
    aspc = [1 (ymax-ymin)/(xmax-xmin) 1];
    % Analytic

    %  Auto-generated by MATLAB on 19-Feb-2021 10:51:09
   figure1 = figure;

    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');

    % Create plot
    plot(xx,yy,'DisplayName','analytical','LineWidth',3,...
        'Color',[0.93 0.69 0.13]);

    plot(x,y,'DisplayName','numerical','MarkerSize',8,'Marker','o',...
        'LineWidth',2,...
        'LineStyle','none',...
        'Color',[0.0 0.45 0.74]);

    % Create ylabel
    ylabel({'$p (Pa)$'},'FontSize',16,'Interpreter','latex');

    % Create xlabel
    xlabel({'$x (m)$'},'FontSize',16,'Interpreter','latex');

    box(axes1,'on');
    % Set the remaining axes properties
    set(axes1,'Color','none','DataAspectRatio',aspc,'FontName',...
        'Times','FontSize',14,'LineWidth',1,'TickLabelInterpreter','latex',...
        'TitleFontWeight','normal','XColor',[0 0 0],'XMinorTick','on','YColor',...
        [0 0 0],'YMinorTick','on','ZColor',[0 0 0]);
    % Uncomment the following line to preserve the X-limits of the axes
    xlim(axes1,[xmin xmax]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim(axes1,[ymin ymax]);
    box(axes1,'on');
    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor','none','FontName','Times','Interpreter','latex');
end
