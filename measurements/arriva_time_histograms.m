clear all
close all



% 10 and 0km fibers
% filename = './20180717-10-0-km-fibers-arrival-times.csv'
% M = csvread(filename);
% figure(1)
% % h = histogram(M(:,2),'BinLimits',[4.90e-05,4.920e-05],'BinWidth',3.449999999970295e-10);
% h = histogram(M(:,2),'BinLimits',[4.90e-05,4.920e-05],'BinWidth',1.2e-9);

% 10 and 10km fibers
% filename = './20180717-10-10-km-fibers-arrival-times.csv'
% M = csvread(filename);
% figure(1)
% % h = histogram(M(:,2),'BinLimits',[-1.00e-7,-7.00e-08,],'BinWidth',3.800000000000338e-11);
%  h = histogram(M(:,2),'BinLimits',[-1.00e-7,-7.00e-08,],'BinWidth',1e-10);

% 5 and 0km fibers
% filename = './20180717-5-0-km-fibers-arrival-times.csv'
% M = csvread(filename);
% figure(1)
% % h = histogram(M(:,2),'BinLimits',[-2.454E-05,-2.450E-05],'BinWidth',5.750000000063430e-11);
% h = histogram(M(:,2),'BinLimits',[-2.454E-05,-2.450E-05],'BinWidth',3e-10);

% 5 and 5km fibers
filename = './20180717-5-5-km-fibers-arrival-times.csv'
M = csvread(filename);
figure(1)
h = histogram(M(:,2),'BinLimits',[-1.6E-08,1.6E-08],'BinWidth',1.310000000000030e-10);
% h = histogram(M(:,2),'BinLimits',[-2.454E-05,-2.450E-05],'BinWidth',3e-10);
set(gca,'YLim',[0 250])





ylabel('Counts');
xlabel('Arrival time (s)');

% Make pretty figure
width = 10;     % Width in inches
height = 6;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 16;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

% print('./10-0-km-fibers.png','-dpng','-r300')
% print('./10-10-km-fibers.png','-dpng','-r300')
% print('./5-0-km-fibers.png','-dpng','-r300')
% print('./5-5-km-fibers.png','-dpng','-r300')