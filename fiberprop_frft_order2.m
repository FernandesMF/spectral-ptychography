clear all

C = 3e+08;  %m/s -> speed of light

L = 10:10:15e+03;     %m
n = 20; %number of samples of the wavefunction

lam0 = 1560e-09;     %m
Dl = 18e-06;    %s/m^2
Dv = -lam0^2/(2*pi*C)*Dl;

sig_wl = [1e-09,3e-09,6e-09,12e-09,15e-09,20e-09,30e-09]; %STDs, in m; std of intensities
sig_f = sig_wl*C/lam0^2;    %Hz
sig_t = 1./(2*pi*sig_f)*sqrt(2);     %s; stds of wavefunctions, not intensities

dt = zeros(length(sig_t),length(L));
a = zeros(length(sig_t),length(L));

for r=1:length(sig_t)
    dt(r,:) = 4*sig_t(r)*sqrt( 1+Dv^2*L.^2/(4*pi^2*sig_t(r)^4) )/n;
    a(r,:) = 2*atan( 4*pi^2*Dv*L./(dt(r,:).^2) )/pi;
end


%% Plot 
wid = 5;    %Width
hei = 3;   %Heigh

nfig = 20;
figure(nfig)
close(nfig)
h = figure(nfig);
hold on
for r = 1:length(sig_t)
    plot(L,a(r,:),'LineWidth',1.5)
end
hold off

hleg = legend(num2str(sig_wl.'),'Location','NorthWest');
htitle = get(hleg,'Title');
set(htitle,'String','Initial wl. band')
title([num2str(n) ' samples inside \pm2\sigma interval'])
xlabel('Fiber Length (m)')
ylabel('Transform order')
set(gca,'FontSize',16)
h.Position = [h.Position(1) h.Position(2) h.Position(1)+100*wid h.Position(2)+100*hei]