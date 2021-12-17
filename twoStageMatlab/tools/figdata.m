function figdata(ref,muY,stdY,tlabel,jump,nf)
    % Create figure
    B=max(max(ref(:,2:end)))*1.2;
    A=min(min(ref(:,2:end)))*0.8;
    C=min(ref(:,1));
    D=max(ref(:,1))*1.01;
    colors = [0.85 0.33 0.10;
        0.07 0.62 1;
        0.93 0.69 0.13;
        0 0 0];

    dasp=[1 1.*(B-A)/(D-C) 200];
    figure1 = figure(nf);
   % Create axes
    axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman',...
        'DataAspectRatio',dasp);
    box(axes1,'on');
    hold(axes1,'all');
    
    for i = 1 : size(muY,2)
        errorbar(ref(1:jump:end,1),muY(1:jump:end,i),stdY(1:jump:end,i),...
            'Parent',axes1,'Color',colors(i,:),'MarkerSize',4,...
            'Marker','o','LineStyle','none','DisplayName',...
            ['mean ' num2str(i,'%d')],'LineWidth',0.5);
    end
    for i = 1 : size(ref,2)-1
        plot(ref(1:end,1),ref(1:end,i+1),'Parent',axes1,'Color',...
            colors(i,:),'MarkerSize',6,'LineWidth',2,...
            'DisplayName',['ref. ' num2str(i,'%d')]);
    end

    xlim(axes1,[0 D])
    ylim(axes1,[A B])
    % Create xlabel
    xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
        'FontAngle','italic','Interpreter','latex');

    % Create ylabel
    ylabel(tlabel,'FontSize',16,'FontName',...
        'Times New Roman','FontAngle','italic','Interpreter','latex');

    % Set the remaining axes properties
    set(axes1,'FontName',...
        'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
        'latex','XMinorTick','on','YMinorTick','on');

    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,'Location','NorthEast','FontSize',8);
    set(legend1,'Box','off');

end

