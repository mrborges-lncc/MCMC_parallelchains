clear all;
[FILENAME, PATHNAME] = uigetfile('../out/mas*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../figuras/';
base_name=[base_aux base_name];
mass = load(line_file);
ma = max(mass);
ma(2) = ma(2)*(1.0+10e-4);
mi = min(mass)*(1.0-10e-4);
A = abs(ma(1)-mi(1));
B = abs(ma(2)-mi(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure1 = figure('XVisual',...
    '0x27 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'FontSize',12,'FontName','Times New Roman',...
    'DataAspectRatio',[1 2*B/A 1]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[mi(1) ma(1)]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[mi(2) ma(2)]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(mass(:,1),mass(:,2),'Parent',axes1,'MarkerSize',8,'Marker','o','LineWidth',2,...
    'Color',[0 0 0]);

% Create xlabel
xlabel('$t$','Interpreter','latex','FontSize',20,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$\theta_{T}$','Interpreter','latex','FontSize',20,...
    'FontName','Times New Roman');

print('-depsc',base_name);
MAX=max(mass(:,2))
MIN=min(mass(:,2))
clear;