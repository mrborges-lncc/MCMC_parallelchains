function plot_rock_poroZinv(rock,g,flag,alpha,rho,titlen,color,lim,vw,n)
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
        'TickLabelInterpreter','latex','XTick',[Lx0:(Lx-Lx0)/4:Lx],...
        'YTick',[Ly0:(Ly-Ly0)/4:Ly],'ZTick',[Lz0:(Lz-Lz0)/5:Lz],...
        'TickLength',[0.0125 0.025],'XMinorTick','on','YMinorTick','on',...
        'ZDir','normal','ZMinorTick','on','YTickLabelRotation',-vw(2),...
        'BoxStyle','full','LineWidth',1,'XTickLabelRotation',vw(1));
    plotCellData(g,perm,'EdgeColor', color,'FaceAlpha', 1.);
    colorbar('horiz','FontSize',12,...
        'TickLabelInterpreter','latex'); axis equal tight; view(vw); box 'on';
    caxis([lim(1) lim(2)]); colormap(jet(55));
    title(titlen,'FontSize',14,'Interpreter','latex');
    % Create labels
    zlabel('$z [mm]$','FontSize',14,'Interpreter','latex');
    ylabel('$y [mm]$','FontSize',14,'Interpreter','latex');
    xlabel('$x [mm]$','FontSize',14,'Interpreter','latex');
    [hc,hh]=colorbarHist(perm,lim,'South', 80);
    pos=get(hc,'Position'); set(hc,'Position',pos - [-.02 -0.03 .05 -.01],'FontSize',11);
    pos=get(hh,'Position'); set(hh,'Position',pos - [-.02 -0.01 .05 -.01]);
    %set(gca,'Position',[.13 .2 .775 .775])
end

