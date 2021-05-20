% load_results.m
% charge les fichiers de sortie.
clear
%clf();

escalax=1;
escalay=1;
filme=10;

disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('../out/b*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
base_aux = '../figuras/mass/';
for i=length(FILENAME):-1:1
    a=FILENAME(i);
    if(a=='.')
        LF=i-1;
    end
end
base_name = [base_aux FILENAME(5:LF)];

time = [];

fid = fopen(line_file,'r');
test=1;
nbstep = 0;
clear satmat
clear velymat
clear velxmat
MIN=1e+28;
MAX=1e-28;
while test==1
    textline = fgetl(fid);
    if textline == -1
        test=0;
    elseif textline(1:10)=='#TIMESTEP '
%    elseif textline(1:10)=='SATURATION'
        nbstep = nbstep+1;
        time = [time str2num(textline(22:length(textline)))];
        mattamp = fscanf(fid,'%f');
	%nbstep
        mattamp = reshape(mattamp,3,length(mattamp)/3)';
        if nbstep==1
            dx = 2.0*mattamp(1,1);
            dy = 2.0*mattamp(1,2);
            nbx = int16((max(mattamp(:,1))/dx)+0.5000000001);
            nby = int16((max(mattamp(:,2))/dy)+0.5000000001);
            xpos = reshape((mattamp(:,1)),nbx,nby,1);
            ypos = reshape((mattamp(:,2)),nbx,nby,1);
       end
        satmat(:,:,nbstep) = reshape(mattamp(:,3),nbx,nby,1);
    elseif textline(1:10)=='VELOCITY -'
        mattamp = fscanf(fid,'%f');
        mattamp = reshape(mattamp,4,length(mattamp)/4)';
        %velxmat(:,nbstep) = mattamp(:,3);
        %velymat(:,nbstep) = mattamp(:,4);
    end
end
fclose(fid);
disp('OK, file loaded')
disp(' ')% satmap.m
% cartes de saturation 3D.

disp('----------------')
disp(' 3D saturation movie ');
MIN=min(min(min((satmat))))
MAX=max(max(max((satmat))))

figure2=figure(...
  'PaperUnits','centimeters',...
  'PaperOrientation','landscape',...
  'PaperSize',[18.92 14.57],...
  'PaperType','<custom>');
%figure(2);
if(filme==1)
    M=moviein(length(time));
end

box on
dx2=dx/2;
dy2=dy/2;
A = (max(max(xpos)))+dx2;
B = (max(max(ypos)))+dy2;
cc= min(min(min(satmat(:,:,:))));
dd= max(max(max(satmat(:,:,:))));
cc=-2.0e-12
dd=6.0e-13
C = cc;
D = 1.0*dd;

for i=48:length(time);
    %cc= (min(min(satmat(:,:,i))));
    %dd= (max(max(satmat(:,:,i))));
    %C = cc;
    %D = 1.01*dd;
    %soma=sum(sum(satmat(:,:,i)))
    for j=1:nby
        for k=1:nbx
            y1 = ypos(k,j)-dy2;
            y2 = ypos(k,j)+dy2;
            x1 = xpos(k,j)-dx2;
            x2 = xpos(k,j)+dx2;
            xx = [x1 x1; x2 x2];
            yy = [y1 y2; y1 y2];
            p = satmat(k,j,i)+[0.0 0.0; 0.0 0.0];
            surf(xx,yy,p);
            hold on
        end
    end
%
    box on
    hold off
    shading('interp');
    alpha(0.85);
%
    axis([0 A 0 B  C D]);
    daspect([1 1 3*abs(D-C)/A]); 
%
    if(dd*cc>=0.0)
        set(gca,'ZTick',[cc 0.5*(dd+cc)  dd],'FontName','times','FontSize',18);
    else
        if((abs(cc)>abs(dd)))
            set(gca,'ZTick',[cc cc/2 0.0  dd],'FontName','times','FontSize',18);
        else
            set(gca,'ZTick',[cc 0.0 dd/2 dd],'FontName','times','FontSize',18);
        end
    end
%
    zlabel ('M_{source} ','FontName','times','FontSize',22);
    set(get(gca,'ZLabel'),'Rotation', 90.0);
    %view(30,15);
    view(55,23);
    %view(0,90);
    %view(140,15);
%
    num=num2str(i-1);
    base_aux = ['_mass-' num];
    time(i)
    base=[base_name base_aux]
%    print('-djpeg90',base)
    print('-depsc','-r100',base)
     if(filme==1)
         M(i)=getframe;
     end
end


if(filme==1)
    movie(M,1,1);
    test=input(' play again (1/0) ');
    while(test==1);
        movie(M,1,3);
        test=input(' play again (1/0) ? ');
    end
    test = input('store it as AVI (1/0) ? ');
    if test==1
%        nomfic = input('file name : ','s');
        nomfic = ['../figuras/movies/' base_name1];
        movie2avi(M,nomfic,'FPS',3);
        movie2avi(M,nomfic,'FPS',4,'quality',100);
    end
end
MIN
MAX
clear
disp(' ')
