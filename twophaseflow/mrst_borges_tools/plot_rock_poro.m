function plot_rock_poro(rock,g,flag,alpha,rho,titlen,color,lim,vw,n)
    fig = figure(n);
    if flag == 'Y'
        perm = (log(rock) - log(alpha))/rho;
    else
        perm = rock;
    end
    plotCellData(g,perm,'EdgeColor', color);
    colorbar('horiz'); axis equal tight; view(vw); box 'on';
    caxis([lim(1) lim(2)]); colormap(jet(55));
    title(titlen,'FontSize',14,'Interpreter','latex');
%     axes('CLim',[0 2],'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',...
%     [50 50 1],'TickDir','both','TickLength',[0.025 0.05],'XMinorTick','on',...
%     'YMinorTick','on','ZDir','reverse','ZMinorTick','on');
    % Create labels
    zlabel('$z (m)$','FontSize',14,'Interpreter','latex');
    ylabel('$y (m)$','FontSize',14,'Interpreter','latex');
    xlabel('$x (m)$','FontSize',14,'Interpreter','latex');
    [hc,hh]=colorbarHist(perm,lim,'South', 80);
    pos=get(hc,'Position'); set(hc,'Position',pos - [.02 0.0 .05 -.01],'FontSize',12);
    pos=get(hh,'Position'); set(hh,'Position',pos - [.02 0.02 .05 -.01]);
    %set(gca,'Position',[.13 .2 .775 .775])
end

