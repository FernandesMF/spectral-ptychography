clear all
close all

% Things to think about: 
%
% fractional FT instead of conventional FT
% initial random estimate: 
%       -non-gaussian complexes
%       -non-smooth
%       -find a haar measure
%       -find that transform to get random polar coordinate numbers
% physical model:
%       -get correct photon wavefunction, right now I have a wf for a massive particle
%       -spectrum filter: does it add phase? Is it just gaussian? do variances increase with tilt?
%       -implement tilt
%       -implement frequency scale (actually, that is easier than I first thought)

% Cool music that played while I was doing this:
% Wild Child - Sinking Ship
% Aquilo - Thin

% PARAMETERS THAT WORKED WELL
% (9jul) T=6,dt=0.1,k0=5,beta=0.7,mu=10:1:50,sig=3,500 ptycho its, PIE update
% (9jul) T=6,T=6,dt=0.1,k0=5,beta=1.0,delta=0.13,mu=10:1:50,sig=3,500 ptycho its, rPIE update

%% Parameters

% Time and frequency vectors vector
% fs = 1e10;
% dt = 1/fs;
% L = 4000;
% t = (-L/2:L/2-1)*dt;
% f = fs*(-L/2:L/2-1)/L;

dt = 0.1;
T = 6;
t = -T/2+dt:dt:T/2;
f = 1:1:length(t);



% Photon gaussian wavefunction
% from https://en.wikipedia.org/wiki/Wave_packet , section about dispersive
% solution; I will use x = 0 further

% lam = 700*1e-09 %(700nm)
% k0 = 2*pi/lam;

k0 = 5;

% Filters
mu = 10:3:50;
sig = ones(size(mu))*0.5e1;
N = length(mu);

% Ptychography/agorithm parameters
a = 0.5;                    %fractional order of the intermediate plane
cycling_method = 'random';  %'random', 'linear' 
stp_crit = 'npty';          %'rel-change', 'npty'
stp_param = 5000;           %either the rel change threshold or the number of ptycho. iterations
npty_max = 5000;
beta = 0.5;                 %step size/correction multiplier
delta = 0.05;               %Wiener filter/rPIE parameter
mom_mode = 'nesterov';       %'none','classic','nesterov'
mom_T = 30;                  %number of iterations to gather when using momentum
eta = 0.9;                  %momentum 'friction' (or better yet, '1-friction'







%% Simulate gaussian wavefunction of a photon
% from https://en.wikipedia.org/wiki/Wave_packet , section about dispersive solution; I will use x=0
% further

psi = 1./sqrt(1+1i*2*t).*exp( -(k0*t).^2./(1+4*t.^2) ).*exp( 1i*(k0^2*t)./(4+16*t.^2) );
psi = ifftshift(psi);
Psi = fft(psi);





%% Apply filters in the frequency domain, generate ptychographic data
% gaussian: g(x) = 1/sig*sqrt(2pi) exp( -(x-mu)^2/(2*sig^2))

% Making gaussian filters
S = zeros(N,length(psi));
for r=1:N
    S(r,:) = exp( -(f-mu(r)).^2/(2*sig(r)^2) )/(sig(r)*sqrt(2*pi));
end
% S = ifftshift(S,2);

% Applying filters to original spectra
Psi_S = zeros(size(S));
for r=1:N
    Psi_S(r,:) = S(r,:).*Psi;
end

% Fourier plane amplitudes
A_F = abs(Psi_S);

% Intermediate plane amplitudes
% fractional FT: https://www.mathworks.com/matlabcentral/fileexchange/41351-frft-m
frPsi_S = zeros(size(Psi_S));
psi_S = zeros(size(Psi_S));
for r=1:N
    psi_S(r,:) = ifft(Psi_S(r,:));
    frPsi_S(r,:) = frft(Psi_S(r,:),-1+a);
end
A_frF = abs(frPsi_S);




%% Run ptychographic iterative engine

% Initial random estimate (smoothed, gaussian complex numbers)
obj = (randn(size(psi))-0.5)+1i*(randn(size(psi))-0.5);
obj_ini = obj;      %saving initial estimate, just for checking purposes...
% obj = smooth(obj,10).';
% obj(200:800) = 0;   %zeroing part of the estimate that should not have any amp
A = sqrt(sum( abs(psi).^2 ));
A_obj = sqrt(sum( abs(obj).^2 ));
obj = A*obj/A_obj;   %normalizing initial estimate
Obj = fft(obj);

% Ptychographic engine
stp_flag = 0;
Ea = [];
Df = [];
Fid = [];
npty = 0;

switch mom_mode
    case 'none'
        disp('Momentum will not be used');
    case {'classic','nesterov'}
        mom_t = 0;
        v = zeros(1,length(Obj));
        Obj_old = Obj;
    otherwise
        error('Unrecognized momentum mode');
end

while(stp_flag==0)
    % Cycling order
    switch cycling_method
        case 'linear'
            s = 1:N;
        case 'random'
            s = randperm(N);
        otherwise
            error('Unrecognized cycling method')
    end
    
    % Ptychographic iteration:
    % -(we wil begin in the Fourier plane)
    % -apply freq. filter
    % -apply inverse frFT (propagate back to intermediate plane)
    % -correct ampitudes (A_frF)
    % -propagate to Fourier plane (apply forward frFT)
    % -correct estimate (Wiener) -> in the Fourier plane... that's where I know the illumination
    for r=1:N
        G = S(s(r),:).*Obj;
        gfr = frft(G,-1+a).';
        gfr_c = A_frF(s(r),:).*exp(1i*angle(gfr));
        G_c = frft(gfr_c,1-a).';
        
        U = updweight(S(s(r),:),delta);
        Obj_new = Obj + beta*U.*( G_c-G );
    end
    
    % Momentum update
    if(strcmp(mom_mode,'classic')||strcmp(mom_mode,'nesterov'))
        mom_t = mom_t+1;
        if(mom_t==mom_T)
            v = eta*v + Obj_new-Obj_old;
            switch mom_mode
                case 'classic'
                    Obj_new =  Obj_old + v;
                case 'nesterov'
                    Obj_new = Obj_new + eta*v;
            end
            mom_t = 0;
            Obj_old = Obj_new;
        end
    end
   
    % Calculate errors, fidelity, rel. change
    ea = sum( (abs(Obj_new)-abs(Psi)).^2 );
    df = sqrt(sum( abs(Obj_new-Obj).^2 ))/sqrt(sum(abs(Obj).^2));
    fid = abs(Obj_new*Psi')/( sqrt(sum(abs(Obj_new).^2))*sqrt(sum(abs(Psi).^2)) );
    
    Ea = [Ea ea];
    Df = [Df df];
    Fid = [Fid fid];
    
    npty = npty+1;
    Obj = Obj_new;
    
    % Check stop condition
    switch stp_crit
        case 'rel-change'
            if(df<stp_param||npty>=npty_max)
                stp_flag = 1;
            end
        case 'npty'
            if(npty>=stp_param)
                stp_flag = 1;
            end
        otherwise
            error('Unrecognized stop criterion')
    end
    
end






%% Plots, show results

% Correct wavefunction (in object and Fourier domain)
figure(1)
subplot(2,2,1)
plot(t,abs(psi))
subplot(2,2,2)
plot(t,angle(psi)/pi)
subplot(2,2,3)
plot(fftshift(abs(Psi)))
subplot(2,2,4)
plot(fftshift(angle(Psi)/pi))

% Frequency filters (and F. domain wf.)
figure(2)
subplot(2,1,1)
plot(fftshift(abs(Psi)))
subplot(2,1,2)
hold on
for r=1:length(mu)
    plot(fftshift(S(r,:),2))
end
hold off

% Fourier domain amplitudes
figure(3)
hold on
for r=1:length(mu)
    plot(fftshift(A_F(r,:)))
end
hold off

% Intermediate plane ampitudes
figure(4)
hold on
for r=1:length(mu)
    plot(A_frF(r,:))
end
hold off

% Initial estimate
figure(5)
subplot(1,2,1)
plot(t,abs(obj_ini),t,abs(obj))
subplot(1,2,2)
plot(t,angle(obj_ini)/pi,t,angle(obj)/pi)

% Final estimate
koo = 1:length(Psi);
figure(6)
subplot(1,2,1)
plot(koo,fftshift(abs(Psi)),'r-',koo,fftshift(abs(Obj_new)),'kx');
subplot(1,2,2)
plot(koo,fftshift(angle(Psi)/pi),'r-',koo,fftshift(angle(Obj_new)/pi),'kx');
figure(7)
subplot(1,2,1)
plot(koo,abs(psi),'r-',koo,abs(ifft(Obj_new)),'kx');
subplot(1,2,2)
plot(koo,angle(psi)/pi,'r-',koo,angle(ifft(Obj_new))/pi,'kx');

% Error curves
figure(8)
subplot(3,1,1)
plot(Ea); ylabel('Ea');
subplot(3,1,2)
plot(Df); ylabel('Df');
subplot(3,1,3)
plot(Fid); ylabel('Fid');

%seminar figures
% figure(10)
% subplot(2,2,1)
% plot(t,ifftshift(abs(psi)),'LineWidth',1.5)
% ylabel('|f|')
% subplot(2,2,2)
% plot(t,ifftshift(angle(psi)/pi),'LineWidth',1.5)
% ylabel('arg(f)')
% subplot(2,2,3)
% plot(fftshift(abs(Psi)),'LineWidth',1.5)
% ylabel('|F|')
% subplot(2,2,4)
% plot(fftshift(angle(Psi)/pi),'LineWidth',1.5)
% ylabel('arg(F)')
% 
% koo = 1:length(Psi);
% figure(11)
% subplot(1,2,1)
% plot(koo,fftshift(abs(Psi)),'r-',koo,fftshift(abs(Obj_new)),'kx','LineWidth',1.5);
% subplot(1,2,2)
% plot(koo,fftshift(angle(Psi)/pi),'r-',koo,fftshift(angle(Obj_new))/pi,'kx','LineWidth',1.5);
% 
% figure(12)
% subplot(2,1,1)
% plot(fftshift(abs(Psi)),'LineWidth',1.5)
% subplot(2,1,2)
% hold on
% for r=1:length(mu)
%     plot(S(r,:),'LineWidth',1.5)
% end
% hold off
% 
% figure(13)
% hold on
% for r=1:length(mu)
%     plot(A_frF(r,:),'LineWidth',1.5)
% end
% hold off