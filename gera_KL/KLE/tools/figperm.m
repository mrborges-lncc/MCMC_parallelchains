function figperm(dout,nx,ny)   
mi=0;
ma=1;
figure1 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',[0.6345 0.6345 11.7 9.7],...
      'PaperSize',[12 10],...
      'PaperType','<custom>');

    axes1 = axes('CLim',[mi ma],...
      'DataAspectRatio',[1 1 2*(ma-mi)/(x*dx)],...%'DataAspectRatio',[x*dx y*dx 3*(ma-mi)/(x*dx)],...
      'FontName','times',...
      'Position',posi,...
      'FontSize',12,...
      'Parent',figure1);
    xlabel(axes1,'$x (m)$','FontName','times','FontSize',16,'Interpreter','latex');
    ylabel(axes1,'$y (m)$','FontName','times','FontSize',16,'Interpreter','latex');
    axis(axes1,[0 x*dx 0 y*dy mi ma]);
    view(axes1,[0 90]);
    %view(axes1,[35 30]);
    box(axes1,'on');
    hold(axes1,'all');
  dy = 1/ny;
  dx = 1/nx;
  figure(1)
  hold on
  for j=ny:-1:1
      yi = (ny-j)*dy;
      yf = yi+dy;
      yy = [yi yf; yi yf];
      for i=1:nx
          xf = i*dx;
          xi = xf - dx;
          xx = [xi xi; xf xf];
          p = dout(j,i)+[0 0; 0 0];
          surf(xx,yy,p);
      end
  end
  view(0,90)

end

