function CreepModeFig(X1, YMatrix1, X2, YMatrix2)
%CREATEFIGURE(X1, YMatrix1, X2, YMatrix2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 01-Dec-2020 12:51:59

% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',subplot1,'LineWidth',2);
set(plot1(1),'DisplayName','Lid: Diffusion','LineStyle','--',...
    'Color',[0 0 1]);
set(plot1(2),'DisplayName','Lid: Basal Slip','LineStyle','--',...
    'Color',[0 1 1]);
set(plot1(3),'DisplayName','Lid: Grain Boundary Sliding','LineStyle','--',...
    'Color',[0 1 0]);
set(plot1(4),'DisplayName','Shell: Diffusion','Color',[0 0 1]);
set(plot1(5),'DisplayName','Shell: Basal Slip','Color',[0 1 1]);
set(plot1(6),'DisplayName','Shell: Grain Boundary Sliding','Color',[0 1 0]);

% Create ylabel
ylabel({'Depth (km)'});

% Create xlabel
xlabel({'Time (kyr)'});

% Create title
title({'Strain Rate: 1x10^1^0 sec^-^1'});

box(subplot1,'on');
axis(subplot1,'ij');
hold(subplot1,'off');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.250964207858093 0.548611195065362 0.191352341601689 0.189114385865271]);

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure1);
hold(subplot2,'on');

% Create multiple lines using matrix input to plot
plot2 = plot(X2,YMatrix2,'Parent',subplot2,'LineWidth',2);
set(plot2(1),'DisplayName','Lid: Diffusion','LineStyle','--',...
    'Color',[0 0 1]);
set(plot2(2),'DisplayName','Lid: Basal Slip','LineStyle','--',...
    'Color',[0 1 1]);
set(plot2(3),'DisplayName','Lid: Grain Boundary Sliding','LineStyle','--',...
    'Color',[0 1 0]);
set(plot2(4),'DisplayName','Shell: Diffusion','Color',[0 0 1]);
set(plot2(5),'DisplayName','Shell: Basal Slip','Color',[0 1 1]);
set(plot2(6),'DisplayName','Shell: Grain Boundary Sliding','Color',[0 1 0]);

% Create ylabel
ylabel({'Depth (km)'});

% Create xlabel
xlabel({'Time (kyr)'});

% Create title
title({'Strain Rate: 3x10^1^0 sec^-^1'});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 500]);
box(subplot2,'on');
axis(subplot2,'ij');
hold(subplot2,'off');
% Create legend
legend2 = legend(subplot2,'show');
set(legend2,...
    'Position',[0.695201021718872 0.696103953233511 0.191352341601687 0.191349921305234]);
