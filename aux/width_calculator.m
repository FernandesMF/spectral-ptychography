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
a = 0.0;            %chirp coeficient (no units)
lam0 = 1560e-09;     %m
sig_wl = 30e-09     %m

% f0 = C/lam0;
f0 = 0;
sig_f = sig_wl*C/lam0^2;    %Hz
sig_t = 1/(2*pi*sig_f)     %s

%remember: psi is sqrt of the prob. dens.
psi = exp(-t.^2/(2*(sig_t*sqrt(2))^2) ).*exp(1i*2*pi*f0.*t).*exp(1i*a*t.^2/(sqrt(2)*sig_t)^2 )/sqrt(sig_t*sqrt(2*pi));
psi = fftshift(psi);
t = fftshift(t);



% Propagation through fiber (does it get wide enough?)
L = 100:100:1.5e+04;     %m
% L = 0;
Dl = 18e-06;    %s/m^2
Dv = -lam0^2/(2*pi*C)*Dl;
v = 1.5e+08;
atten = 1;

sig_t2 = sqrt( sum(t.^2.*abs(psi).^2*dt) - sum(t.*abs(psi).^2*dt)^2 )

SigsL_t = zeros(1,length(L));
for r=1:length(L)
    [psiL,H] = fiberprop(lam0,t,psi,L(r),Dv,v,atten);
    sigL_t = sqrt( sum(t.^2.*abs(psiL).^2*dt) - sum(t.*abs(psiL).^2*dt)^2 );
    SigsL_t(r) = sigL_t;    %SigsL_t has stds of INTENSITY PROFILE
end

sprintf('%g\t%g\n',reshape([L; SigsL_t],[1 2*length(L)]))   % hell yeah







%% Plots, show results

figure(1)
plot(L,SigsL_t,'-o');


L = 100:100:1.5e+04;
