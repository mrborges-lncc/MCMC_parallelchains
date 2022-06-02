function figerrorColor2(X1, Y1, Y2, Y3, flag, maxi, maxi2)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data
fac = 1;
if(flag==1)
    %  Auto-generated by MATLAB on 19-Dec-2011 08:09:49
    B=max(Y1)*1.1;
    A=0.0;%min(Y1);
    if(maxi>0)
        B=maxi;
    end
    C=min(X1);
    D=max(X1);

    % Create figure
    figure1 = figure(56);

    % Create axes
    dasp=[1 (B-A)/(D-C) 200];
%     axes1 = axes('Parent',figure1,'YScale','linear','YMinorTick','on',...
    axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
        'XMinorTick','on','FontSize',14,'TickDir','both',...
        'FontName','Times New Roman','TickLabelInterpreter','latex',...
        'DataAspectRatio',dasp);
    box(axes1,'on');
    hold(axes1,'all');
    xlim(axes1,[0 D])
    ylim(axes1,[A B])

    % Create semilogy
    % plot(X1,Y1,'Parent',axes1,'MarkerSize',9,'Marker','o','LineWidth',2,...
    %     'DisplayName','Independent sampling',...
    %     'Color',[0 0 0]);
    semilogy(X1,Y1,'Parent',axes1,'MarkerSize',9,'Marker','none','LineWidth',2,...
        'DisplayName','Independent sampling',...
        'Color',[0 0 0]);

    % Create xlabel
    xlabel('Accepted iterations','FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman');

    % Create ylabel
    % ylabel('Fractional flow error','FontWeight','bold','FontSize',16,...
    %     'FontName','Times New Roman');
    ylabel('Total relative error','FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman');

    % Create legend
    % legend(axes1,'show');
end
if(flag==0)
    %  Auto-generated by MATLAB on 19-Dec-2011 08:09:49
    B=max(maxi2,maxi);%max(max(Y2),max(Y3))*1.1;
    A=min(min(Y3),min(min(Y2),min(Y2)));
%     if(maxi>0)
%         B=maxi;
%     end
    C=min(X1);
    D=max(X1);

    % Create figure
    figure1 = figure(456);

    % Create axes
    dasp=[fac (B-A)/(D-C) 200];
%     axes1 = axes('Parent',figure1,'YScale','linear',...
    axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
        'TickDir','both','Position',[0.11 0.13 0.8 0.815],...
        'YColor',[0 0 1],'LineWidth',1,...
        'YMinorTick','on','TickLabelInterpreter','latex',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'DataAspectRatio',dasp);
    box(axes1,'on');
    hold(axes1,'all');
    xlim(axes1,[C D])
    ylim(axes1,[A B])

    % Create semilogy
    % plot(X1,Y1,'Parent',axes1,'MarkerSize',9,'Marker','o','LineWidth',2,...
    %     'DisplayName','Independent sampling',...
    %     'Color',[0 0 0]);
    plot(X1,Y2,'Parent',axes1,'MarkerSize',9,'Marker','none',...
        'LineWidth',2,'LineStyle','-',...
        'Color',[0 0 1]);

    % Create xlabel
    xlabel('$\mathsf{MCMC}$ iterations','FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman','Interpreter','latex');

    % Create ylabel
    % ylabel('Fractional flow error','FontWeight','bold','FontSize',16,...
    %     'FontName','Times New Roman');
    ylabel('Relative error ($p$)','Interpreter','latex',...
        'FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman','Color',[0 0 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     B=maxi;
    dasp2=[fac (B-A)/(D-C) 200];
    axes2 = axes('Parent',figure1,'YScale','log',...
        'YAxisLocation','right','XColor','none',...
        'Position',[0.11 0.13 0.8 0.815],...
        'YColor',[1 0 0],'LineWidth',1.5,...
        'YMinorTick','on',...
        'YMinorTick','on','TickLabelInterpreter','latex',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'DataAspectRatio',dasp2,'Color','none');
%     box(axes2,'on');
    hold(axes2,'all');
    xlim(axes2,[C D])
    ylim(axes2,[A B])

    plot(X1,Y3,'Parent',axes2,'MarkerSize',9,'Marker','none',...
        'LineWidth',2,'LineStyle','-',...
        'Color',[1 0 0]);
    ylabel('Relative error ($\mathsf{F}$)','Interpreter','latex',...
        'FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman');

    % Create legend
%     legend1=legend(axes1,'show');
%     set(legend1,'FontAngle','italic');
%     set(legend1,'box','off');
end
 if(flag==2)
    %  Auto-generated by MATLAB on 19-Dec-2011 08:09:49
    B=maxi2;%max(max(Y2),max(Y3))*1.1;
    A=0.0;%min(Y1);
%     if(maxi>0)
%         B=maxi;
%     end
    C=min(X1);
    D=max(X1);

    % Create figure
    figure1 = figure();

    % Create axes
    dasp=[1 (B-A)/(D-C) 200];
    %axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    axes1 = axes('Parent',figure1,'YScale','linear',...
    'Position',[0.11 0.13 0.8 0.815],...
        'YColor',[0 0 1],...
        'YMinorTick','on',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'DataAspectRatio',dasp);
    box(axes1,'on');
    hold(axes1,'all');
    xlim(axes1,[0 D])
    ylim(axes1,[A B])

    % Create semilogy
    % plot(X1,Y1,'Parent',axes1,'MarkerSize',9,'Marker','o','LineWidth',2,...
    %     'DisplayName','Independent sampling',...
    %     'Color',[0 0 0]);
    plot(X1,Y2,'Parent',axes1,'MarkerSize',9,'Marker','none',...
        'LineWidth',2,'LineStyle','-',...
        'Color',[0 0 1]);

    % Create xlabel
    xlabel('$\mathsf{MCMC}$ iterations','FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman','Interpreter','latex');

    % Create ylabel
    % ylabel('Fractional flow error','FontWeight','bold','FontSize',16,...
    %     'FontName','Times New Roman');
    ylabel('Relative error ($\mathsf{F}$)','Interpreter','latex',...
        'FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman','Color',[0 0 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B=maxi;
    dasp2=[1 (B-A)/(D-C) 200];
    axes2 = axes('Parent',figure1,'YScale','linear',...
        'YAxisLocation','right',...
        'Position',[0.11 0.13 0.8 0.815],...
        'YColor',[1 0 0],...
        'YMinorTick','on',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'DataAspectRatio',dasp2,'Color','none');
%     box(axes2,'on');
    hold(axes2,'all');
    xlim(axes2,[0 D])
    ylim(axes2,[A B])

    plot(X1,Y3,'Parent',axes2,'MarkerSize',9,'Marker','none',...
        'LineWidth',2,'LineStyle','-',...
        'Color',[1 0 0]);
%     ylabel('Relative error (${u_y}$)','Interpreter','latex',...
%         'FontWeight','bold','FontSize',16,...
%         'FontName','Times New Roman');

    % Create legend
%     legend1=legend(axes1,'show');
%     set(legend1,'FontAngle','italic');
%     set(legend1,'box','off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B=maxi;
    dasp2=[1 (B-A)/(D-C) 200];
    axes2 = axes('Parent',figure1,'YScale','linear',...
        'YAxisLocation','right',...
        'Position',[0.11 0.13 0.8 0.815],...
        'YColor',[1 0 0],...
        'YMinorTick','on',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'DataAspectRatio',dasp2,'Color','none');
%     box(axes2,'on');
    hold(axes2,'all');
    xlim(axes2,[0 D])
    ylim(axes2,[A B])

    plot(X1,Y1,'Parent',axes2,'MarkerSize',9,'Marker','none',...
        'LineWidth',2,'LineStyle','--',...
        'Color',[1 0 0]);
    ylabel('Relative error ($u_y,\ {p}$)','Interpreter','latex',...
        'FontWeight','bold','FontSize',16,...
        'FontName','Times New Roman');

    % Create legend
%     legend1=legend(axes1,'show');
%     set(legend1,'FontAngle','italic');
%     set(legend1,'box','off');
end
   
