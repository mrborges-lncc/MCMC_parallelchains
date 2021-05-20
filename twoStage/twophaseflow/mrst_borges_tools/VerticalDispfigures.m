function np = VerticalDispfigures(uu,G,W,printa,vw,nome,et,num,nprj,t,npk,ndt,lim)
    t = ['$\quad t = ' num2str(t,'%5.3f') '\ days$'];
    Lz0 = min(G.nodes.coords(:,3));
    Lz  = max(G.nodes.coords(:,3));
    Z   = Lz-Lz0;
    if printa == 1 && mod(num-ndt,nprj) == 0 && num > ndt
        figure(6);
        plotNodeDataDeformed(G, uu(:,3), uu, ...
            'FaceAlpha', 0.95, 'EdgeAlpha', 0.4, 'EdgeColor', 'none');
        plotWell(G, W, 'height', Z*0.2,'FontSize',10,'Interpreter','latex')
        [hc] = colorbar('horiz'); colormap(jet(55)); axis equal tight; view(vw);
        set(gca,'ZDir', 'reverse'), title(['Vertical disp. [$m$]' t],...
            'FontSize',12,'Interpreter','latex',...
            'Position', [20 5.5 1150]);
        set(gca, 'DataAspect',[1 1 1],'Color','none','Box','on')
        caxis(lim)
        % Create labels
        zlabel('$z (m)$','FontSize',14,'Interpreter','latex');
        ylabel('$y (m)$','FontSize',14,'Interpreter','latex');
        xlabel('$x (m)$','FontSize',14,'Interpreter','latex');
        base=['figuras/zdisp_' nome '-' num2str(npk,'%d')];
        set(gcf,'PaperPositionMode','auto');
        print('-dpng','-r600', base);
        pause(et);
       clf; close all;
        figure(7);
        plotGridDeformed(G, uu*10,...
            'FaceAlpha', 0.95, 'EdgeAlpha', 0.4, 'EdgeColor', 'k');
        plotWell(G, W, 'height', Z*0.2,'FontSize',10,'Interpreter','latex')
        %colorbar('horiz'); 
        colormap(jet(55)); axis equal tight; view(vw);
        set(gca,'ZDir', 'reverse'), title(['Deformed grid, ' t],...
            'FontSize',12,'Interpreter','latex',...
            'Position', [19 6.7 1150]);

        set(gca, 'DataAspect',[1 1 1/3],'Color','none','Box','on')
        % Create labels
        zlabel('$z (m)$','FontSize',14,'Interpreter','latex');
        ylabel('$y (m)$','FontSize',14,'Interpreter','latex');
        xlabel('$x (m)$','FontSize',14,'Interpreter','latex');
        %zlim([(Lz0-(Lz-Lz0)*0.00) Lz]);
        base=['figuras/deformz_' nome '-' num2str(npk,'%d')];
        set(gcf,'PaperPositionMode','auto');
        print('-dpng','-r600', base);
        pause(et);
        clf; close all;
        np = npk + 1;
    else
        np = npk;
    end
end

