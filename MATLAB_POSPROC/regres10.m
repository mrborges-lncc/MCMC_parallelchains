function [fitResults1] =regres10(X1, Y1)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 03-Dec-2018 17:47:58

% Create figure
% figure1 = figure;

% Create axes
% axes1 = axes('Parent',figure1);
% box(axes1,'on');
% hold(axes1,'all');

% Create plot
% plot1 = plot(X1,Y1,'Parent',axes1,'Marker','o','LineStyle','none',...
%     'DisplayName','data 1');

% Get xdata from plot
xdata1 = X1;%get(plot1, 'xdata');
% Get ydata from plot
ydata1 = Y1;%get(plot1, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);

% Remove NaN values and warn
nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask1) = [];
    ydata1(nanMask1) = [];
end

% Find x values for plotting the fit based on xlim
% axesLimits1 = xlim(axes1);
% xplot1 = linspace(axesLimits1(1), axesLimits1(2));

% Preallocate for "Show equations" coefficients
coeffs1 = cell(1,1);


fitResults1 = polyfit(xdata1, ydata1, 10);
% Evaluate polynomial
% yplot1 = polyval(fitResults1, xplot1);

% Save type of fit for "Show equations"
fittypesArray1(1) = 11;

% Save coefficients for "Show Equation"
coeffs1{1} = fitResults1;

% Plot the fit
% fitLine1 = plot(xplot1,yplot1,'DisplayName','   10th degree','Parent',axes1,...
%     'Tag','10th degree',...
%     'Color',[0.75 0 0.75]);

% Set new line in proper position
% setLineOrder(axes1, fitLine1, plot1);

% "Show equations" was selected
% showEquations(fittypesArray1, coeffs1, 2, axes1);

% Create legend
% legend(axes1,'show');

%-------------------------------------------------------------------------%
function setLineOrder(axesh1, newLine1, associatedLine1)
%SETLINEORDER(AXESH1,NEWLINE1,ASSOCIATEDLINE1)
%  Set line order
%  AXESH1:  axes
%  NEWLINE1:  new line
%  ASSOCIATEDLINE1:  associated line

% Get the axes children
hChildren = get(axesh1,'Children');
% Remove the new line
hChildren(hChildren==newLine1) = [];
% Get the index to the associatedLine
lineIndex = find(hChildren==associatedLine1);
% Reorder lines so the new line appears with associated data
hNewChildren = [hChildren(1:lineIndex-1);newLine1;hChildren(lineIndex:end)];
% Set the children:
set(axesh1,'Children',hNewChildren);

%-------------------------------------------------------------------------%
function showEquations(fittypes1, coeffs1, digits1, axesh1)
%SHOWEQUATIONS(FITTYPES1,COEFFS1,DIGITS1,AXESH1)
%  Show equations
%  FITTYPES1:  types of fits
%  COEFFS1:  coefficients
%  DIGITS1:  number of significant digits
%  AXESH1:  axes

n = length(fittypes1);
txt = cell(length(n + 1) ,1);
txt{1,:} = ' ';
for i = 1:n
    txt{i + 1,:} = getEquationString(fittypes1(i),coeffs1{i},digits1,axesh1);
end
text(.05,.95,txt,'parent',axesh1, ...
    'verticalalignment','top','units','normalized');

%-------------------------------------------------------------------------%
function [s1] = getEquationString(fittype1, coeffs1, digits1, axesh1)
%GETEQUATIONSTRING(FITTYPE1,COEFFS1,DIGITS1,AXESH1)
%  Get show equation string
%  FITTYPE1:  type of fit
%  COEFFS1:  coefficients
%  DIGITS1:  number of significant digits
%  AXESH1:  axes

if isequal(fittype1, 0)
    s1 = 'Cubic spline interpolant';
elseif isequal(fittype1, 1)
    s1 = 'Shape-preserving interpolant';
else
    op = '+-';
    format1 = ['%s %0.',num2str(digits1),'g*x^{%s} %s'];
    format2 = ['%s %0.',num2str(digits1),'g'];
    xl = get(axesh1, 'xlim');
    fit =  fittype1 - 1;
    s1 = sprintf('y =');
    th = text(xl*[.95;.05],1,s1,'parent',axesh1, 'vis','off');
    if abs(coeffs1(1) < 0)
        s1 = [s1 ' -'];
    end
    for i = 1:fit
        sl = length(s1);
        if ~isequal(coeffs1(i),0) % if exactly zero, skip it
            s1 = sprintf(format1,s1,abs(coeffs1(i)),num2str(fit+1-i), op((coeffs1(i+1)<0)+1));
        end
        if (i==fit) && ~isequal(coeffs1(i),0)
            s1(end-5:end-2) = []; % change x^1 to x.
        end
        set(th,'string',s1);
        et = get(th,'extent');
        if et(1)+et(3) > xl(2)
            s1 = [s1(1:sl) sprintf('\n     ') s1(sl+1:end)];
        end
    end
    if ~isequal(coeffs1(fit+1),0)
        sl = length(s1);
        s1 = sprintf(format2,s1,abs(coeffs1(fit+1)));
        set(th,'string',s1);
        et = get(th,'extent');
        if et(1)+et(3) > xl(2)
            s1 = [s1(1:sl) sprintf('\n     ') s1(sl+1:end)];
        end
    end
    delete(th);
    % Delete last "+"
    if isequal(s1(end),'+')
        s1(end-1:end) = []; % There is always a space before the +.
    end
    if length(s1) == 3
        s1 = sprintf(format2,s1,0);
    end
end

