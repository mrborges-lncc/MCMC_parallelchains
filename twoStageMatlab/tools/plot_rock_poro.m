function plot_rock_poro(rock,g,flag,alpha,rho,titlen,color,lim,vw,n)
    fig = figure(n);
    TOL = 1e-8;
    if flag == 'Y'
        perm = (log(rock) - log(alpha))/rho;
    else
        perm = rock;
    end
    if abs(lim(1) - lim(2)) < TOL
        lim = [min(perm) max(perm)];
        if abs(lim(1) - lim(2)) < TOL
            aux = abs(max(perm))*0.1;
            lim = [(min(perm)-aux) (max(perm)+aux)];
        end
    end
    Lx0 = min(g.nodes.coords(:,1));
    Lx  = max(g.nodes.coords(:,1));
    Ly0 = min(g.nodes.coords(:,2));
    Ly  = max(g.nodes.coords(:,2));
    Lz0 = min(g.nodes.coords(:,3));
    Lz  = max(g.nodes.coords(:,3));

    axes('DataAspectRatio',[1 1 1],'FontName','Times','FontSize',12,...
        'PlotBoxAspectRatio',g.cartDims,'TickDir','both',...
        'TickLabelInterpreter','latex','XTick',[Lx0:100:Lx],...
        'YTick',[Ly0:100:Ly],'ZDir','reverse','ZTick',[Lz0 Lz],...
        'TickLength',[0.0125 0.025],'XMinorTick','on','YMinorTick','on',...
        'ZDir','reverse','ZMinorTick','on');
    plotCellData(g,perm,'EdgeColor', color);
    colorbar('horiz'); axis equal tight; view(vw); box 'on';
    caxis([lim(1) lim(2)]); colormap(jet(55));
    title(titlen,'FontSize',14,'Interpreter','latex');
    % Create labels
    zlabel('$z (m)$','FontSize',14,'Interpreter','latex');
    ylabel('$y (m)$','FontSize',14,'Interpreter','latex');
    xlabel('$x (m)$','FontSize',14,'Interpreter','latex');
    [hc,hh]=colorbarHist(perm,lim,'South', 80);
    pos=get(hc,'Position'); set(hc,'Position',pos - [-.02 -0.10 .05 -.01],'FontSize',12);
    pos=get(hh,'Position'); set(hh,'Position',pos - [-.02 -0.08 .05 -.01]);
    %set(gca,'Position',[.13 .2 .775 .775])
end

