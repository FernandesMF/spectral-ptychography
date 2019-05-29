%L   1nm        3nm         6nm	        12nm         30nm

C = 3e+08;  %m/s -> speed of light
lam0 = 1560e-09;     %m
Dl = 18e-06;    %s/m^2
Dv = -lam0^2/(2*pi*C)*Dl;

L = 0:10:15e03; %m
sig1_wl = [1e-09,3e-09,6e-09,12e-09,30e-09];
sig1_f = sig1_wl*C/lam0^2;    %Hz
sig1_t = 1./(2*pi*sig1_f)*sqrt(2)     %s; sig1 contains stds of WAVEFUNCTION, not measured intensities
Data= zeros(length(L),6);
Data(:,1) = L;

for r=1:5
    Data(:,r+1) = sig1_t(r)*sqrt( 1 + Dv^2*L.^2/(4*pi^2*sig1_t(r)^4) )/sqrt(2); %Data cointains stds of INTENSITIES
end


figure(21)
plot(Data(:,1),Data(:,2),...
     Data(:,1),Data(:,3),...
     Data(:,1),Data(:,4),...
     Data(:,1),Data(:,5),...
     Data(:,1),Data(:,6))
 legend('1nm','3nm','6nm','12nm','30nm','Location','NorthWest')

 

%% 25ps jitter
 
foo50 = ones(size(Data(:,1)))*50e-12;
foo35 = ones(size(Data(:,1)))*35e-12;
[~,find50] = min(abs(Data(:,2:6)-50e-12));
[~,find35] = min(abs(Data(:,2:6)-35e-12));
P50 = []; P35 = [];
for r=1:length(find50)
    P50 = [P50; Data(find50(r),1),Data(find50(r),r+1)];
    P35 = [P35; Data(find35(r),1),Data(find35(r),r+1)];
end
 
figure(22)
plot(Data(:,1),Data(:,2),...
     Data(:,1),Data(:,3),...
     Data(:,1),Data(:,4),...
     Data(:,1),Data(:,5),...
     Data(:,1),Data(:,6),...
     Data(:,1),foo50,'k--',...
     Data(:,1),foo35,'k--',...
     P50(:,1),P50(:,2),'ko',...
     P35(:,1),P35(:,2),'ko','MarkerSize',10,'Linewidth',1.5);
set(gca,'YLim',[0 110e-12])
legend('1nm','3nm','6nm','12nm','30nm','Location','SouthEast')
xlabel('Propagation length (m)')
ylabel('Pulse width (s)')
title('25ps jitter')

axPos = get(gca,'Position');
xMinMax = xlim;
yMinMax = ylim;
for r=1:length(find50)
    xPlot = P50(r,1); yPlot = P50(r,2);
    xAnnotation = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAnnotation = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
    dim = [xAnnotation yAnnotation .03 .03];
    str = [num2str(P50(r,1)) 'm'];
    if(xAnnotation>=0 && xAnnotation<=1 && yAnnotation>=0 && yAnnotation<=1)
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end
for r=1:length(find35)
    xPlot = P35(r,1); yPlot = P35(r,2);
    xAnnotation = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAnnotation = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
    dim = [xAnnotation yAnnotation .03 .03];
    str = [num2str(P35(r,1)) 'm'];
    if(xAnnotation>=0 && xAnnotation<=1 && yAnnotation>=0 && yAnnotation<=1)
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end

% Printing pretty fig
figure(22)
width = 12;     % Width in inches
height = 6;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 2.5;      % LineWidth
msz = 10;       % MarkerSize'pos = get(gcf, 'Position');
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

% print('25ps-jitter2.png','-dpng','-r300')



%% 5ps jitter

foo10 = ones(size(Data(:,1)))*10e-12;
foo07 = ones(size(Data(:,1)))*07e-12;
[~,find10] = min(abs(Data(:,2:6)-10e-12));
[~,find07] = min(abs(Data(:,2:6)-07e-12));
P10 = []; P07 = [];
for r=1:length(find50)
    P10 = [P10; Data(find10(r),1),Data(find10(r),r+1)];
    P07 = [P07; Data(find07(r),1),Data(find07(r),r+1)];
end

figure(23)
plot(Data(:,1),Data(:,2),...
     Data(:,1),Data(:,3),...
     Data(:,1),Data(:,4),...
     Data(:,1),Data(:,5),...
     Data(:,1),Data(:,6),...
     Data(:,1),foo10,'k--',...
     Data(:,1),foo07,'k--',...
     P10(:,1),P10(:,2),'ko',...
     P07(:,1),P07(:,2),'ko','MarkerSize',10,'Linewidth',1.5);
set(gca,'YLim',[0 15e-12])
set(gca,'XLim',[0 5000])
legend('1nm','3nm','6nm','12nm','30nm','Location','SouthEast')
xlabel('Propagation length (m)')
ylabel('Pulse width (s)')
title('5ps jitter')

axPos = get(gca,'Position');
xMinMax = xlim;
yMinMax = ylim;
for r=1:length(find10)
    xPlot = P10(r,1); yPlot = P10(r,2);
    xAnnotation = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAnnotation = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
    dim = [xAnnotation yAnnotation .03 .03];
    str = [num2str(P10(r,1)) 'm'];
    if(xAnnotation>=0 && xAnnotation<=1 && yAnnotation>=0 && yAnnotation<=1)
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end
for r=1:length(find07)
    xPlot = P07(r,1); yPlot = P07(r,2);
    xAnnotation = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAnnotation = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
    dim = [xAnnotation yAnnotation .03 .03];
    str = [num2str(P07(r,1)) 'm'];
    if(xAnnotation>=0 && xAnnotation<=1 && yAnnotation>=0 && yAnnotation<=1)
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end

figure(23)
width = 12;     % Width in inches
height = 6;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 5.0;      % LineWidth
msz = 10;       % MarkerSize'pos = get(gcf, 'Position');
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

% print('5ps-jitter2.png','-dpng','-r300')