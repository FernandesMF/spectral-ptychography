C = 3e+08;  %m/s -> speed of light

L = 0:10:10e+03;     %m

lam0 = 1560e-09;     %m
Dl = 18e-06;    %s/m^2
Dv = -lam0^2/(2*pi*C)*Dl;
s = [3e-12,5e-12,7e-12,1e-11,1.4e-11,2e-11];

for r=1:length(s)
    af(r,:) = 2*atan(4*pi^2*Dv*L/s(r)^2)/pi;
end

nfig = 10;
figure(nfig)
close(nfig);
figure(nfig)
hold on
for r=1:length(s)
    plot(L,af(r,:),'LineWidth',1.5)
end
xlabel('Prop. distance (m)')
ylabel('Fractional order')
hleg = legend(num2str(s'));
htitle = get(hleg,'Title');
set(htitle,'String','Time resolution (s)')
hold off

%% Pretty fig
% h = figure(nfig)
