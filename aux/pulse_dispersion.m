clear all
% close all

C = 3e+08;  %m/s, speed of light

%% Parameters

% Time and frequency vectors 

dt = 1e-15;        %s
T = 2e-09;          %s
t = -T/2:dt:T/2;    %s

N = length(t);

Fs = 1/dt;          %Hz
f = (0:N-1)*Fs/N;     %Hz



% Photon gaussian wavefunction
% gaussian with quadratic phase profile (which can be pretty much any... this is controllable)
a = 0.001;            %chirp coeficient (no units)
lam0 = 1560e-09;     %m
sig_wl = 3e-09     %m

% f0 = C/lam0;
f0 = 0;
sig_f = sig_wl*C/lam0^2;    %Hz
sig_t = 1/(2*pi*sig_f)     %s

%remember: psi is sqrt of the prob. dens.
psi = exp(-t.^2/(2*(sig_t*sqrt(2))^2) ).*exp(1i*2*pi*f0.*t).*exp(1i*a*t.^2/(sqrt(2)*sig_t)^2 )/sqrt(sig_t*sqrt(2*pi));
psi = fftshift(psi);
t = fftshift(t);



% Propagation through fiber (does it get wide enough?)
L = 5e+02;     %m
% L = 0;
Dl = 18e-06;    %s/m^2
Dv = -lam0^2/(2*pi*C)*Dl;
v = 1.5e+08;
atten = 1;
[psiL,H] = fiberprop(lam0,t,psi,L,Dv,v,atten);

sig_t2 = sqrt( sum(t.^2.*abs(psi).^2*dt) - sum(t.*abs(psi).^2*dt)^2 )
sigL_t2 = sqrt( sum(t.^2.*abs(psiL).^2*dt) - sum(t.*abs(psiL).^2*dt)^2 )
% sig_tL = std(psiL)



%% Plots, show results

figure(1)
subplot(2,2,1)
plot(t,abs(psi));
subplot(2,2,2)
plot(t,angle(psi))
subplot(2,2,3)
plot(t,abs(psiL));
subplot(2,2,4)
plot(t,angle(psiL))

b = Dv*L/pi;
figure(2)
subplot(1,2,1)
plot(angle(H))
subplot(1,2,2)
plot(abs(fft(psi)));