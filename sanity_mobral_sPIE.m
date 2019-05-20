% This routine will just test if the spectralPIE routine behaves accordingly when the
% input is actually just the correct function. sPIE should just leave the estimate as
% it is, since there will be no correction. Right?

% I guess this is just a sanity test based on realistic_ptycho_sim routine, but without
% resampling.

clear all
close all

%% Parameters

% ATTENTION: use (nano) SI units throughout this code
C = 3e+08;  % nm/ns, speed of light

% Mock time/ frequency domains
N	= 10;
t   = 1:N;
f   = 1:N;

% Mock interference filter
par_if.nthetas	= N;             % in radians
par_if.r        = 3;

% Ptychography algorithm parameters
par_ptycho.init_est_method  = 'random';
par_ptycho.cycling_method   = 'random';  %'random', 'linear' 
par_ptycho.stp_crit     = 'npty';        %'rel-change', 'npty'
par_ptycho.dthresh      = 1e-02;         % either the rel change threshold or the number of ptycho. iterations
par_ptycho.npty_max     = 500;
par_ptycho.n_rerand     = 0;          
par_ptycho.beta     = 0.5;               % step size/correction multiplier
par_ptycho.delta    = 0.05;              % Wiener filter/rPIE parameter
% par_ptycho.Nmasks   = length(par_if.thetas);    % number of spectral masks to be used
par_ptycho.Nmasks   = par_if.nthetas;

par_ptycho.mom_mode     = 'none';    %'none','classic','nesterov'
par_ptycho.mom_T    = 30;                % number of iterations to gather when using momentum
par_ptycho.eta      = 0.9;               % momentum 'friction' (or better yet, '1-friction')

% par_ptycho.frft_a   = 2*atan(4*pi^2*par_prop.Dv*par_prop.L/(par_detec.timeres^2))/pi; %FIX: change dt to detection dt
par_ptycho.frft_a   = 0.3;

% Photon gaussian wavefunction
% gaussian with quadratic phase profile (which can be pretty much any... this is controllable)
par_wf.a   = 0.1;               % chirp coeficient (no units)
par_wf.lam0    = 1560e-00;      % nm
par_wf.sig_wl  = 3e-00;         % nm

par_wf.f0  = C/par_wf.lam0;     % GHz
par_wf.sig_f   = par_wf.sig_wl*C/par_wf.lam0^2;     % GHz
par_wf.sig_t   = 1/(2*pi*par_wf.sig_f);             % ns

% Propagation through fiber
par_prop.L  = 100;     % m
par_prop.Dl = 18e-06;   % s/m^2
par_prop.Dv = -par_wf.lam0^2/(2*pi*C)*par_prop.Dl;
par_prop.v  = 1.5e+08;


%% Generate photon initial state

psi     = rand(1,N)-0.5 + 1i*(rand(1,N)-0.5);
Psi     = fft(psi);



%% Generate spectral masks
I = zeros(length(par_if.nthetas),N);
% I = ones(size(I)); % for debugging purposes...
I(1,1:par_if.r)     = 1;
for s=2:par_if.nthetas
    I(s,:)  = circshift(I(s-1,:),[0 1]);
end



%% Propagate filtered wavefunctions
Psi_filt    = zeros(par_if.nthetas,N);
Psi_L       = zeros(par_if.nthetas,N);
psi_L       = zeros(par_if.nthetas,N);
for s=1:par_if.nthetas
    Psi_filt(s,:)   = I(s,:).*Psi;
%     %[Psi_L(s,:),H]  = fiberprop(par_wf.lam0,t,Psi_filt(s,:),par_prop.L,par_prop.Dv,par_prop.v,1);
%     
    %f0  = C/par_wf.lamb0;
    f0  = 0;
    td  = par_prop.L/par_prop.v;
    b   = par_prop.Dv*par_prop.L/pi;
    %H   = atten*exp(-1i*( 2*pi*td*(f-f0) + pi^2*b*(f-f0) ));
    H   = exp(-1i*( pi^2*b*(f-f0).^2 ));
    Psi_L(s,:)  = Psi_filt(s,:).*H;
    psi_L(s,:)  = ifft(Psi_L(s,:));

%     psi_L(s,:)  = frft(Psi_filt(s,:),-1+par_ptycho.frft_a);
%     Psi_L(s,:)  = frft(Psi_filt(s,:),par_ptycho.frft_a);
%     psi_L(s,:)  = ifft(Psi_L(s,:));
end



%% Run ptychography
%FIX: implement rerandomization
disp('Starting Ptychography')

% Initial random estimate (smoothed, gaussian complex numbers)
% obj     = (randn(1,rs_N)-0.5)+1i*(randn(1,rs_N)-0.5);
obj = psi;           % sanity test
obj_ini	= obj;          % saving initial estimate, just for checking purposes...
% obj = smooth(obj,10).';


A	= sqrt(sum( abs(psi).^2 ));
A_obj   = sqrt(sum( abs(obj).^2 ));
obj = A*obj/A_obj;      % normalizing initial estimate
% Obj = fft(obj);
Obj = Psi;           % sanity test
Obj_ini = Obj;          % saving initial estimate for checking purposes

% Ptychographic engine
% [Obj,Ea,Df,Fid,bad_res] =  spectralPIE(par_ptycho,Obj,I,abs(psi_L),Psi);
[Obj,Ea,Df,Fid,bad_res] =  spectralPIE_prop(par_ptycho,par_prop,f,Obj,I,abs(psi_L),Psi);





%% Show results

% Original wavefunctions
figure(1)
subplot(2,1,1); plot(t,abs(psi)); xlabel('t'); ylabel('|\psi(t)|');
subplot(2,1,2); plot(f,abs(Psi)); xlabel('f'); ylabel('|\Psi(f)|');

% Filter transmittances and filtered wavefunction (freq. domain)
figure(2)
subplot(3,1,1); plot(f,abs(Psi)); xlabel('f'); ylabel('|\Psi(f)|');
subplot(3,1,2); plot(f,abs(I)); xlabel('f'); ylabel('|T_e(f)|');
% plot(f,abs(I),f,angle(I))
subplot(3,1,3); plot(f,abs(Psi_filt)); xlabel('f'); ylabel('|\Psi_{\rm filt}|');

% % Original and filtered wavefunctions (time and freq. domains)
figure(3)
subplot(3,1,1); plot(t,abs(psi)); xlabel('t'); ylabel('|\psi(t)|');
subplot(3,1,2); plot(t,abs(psi_L)); xlabel('t'); ylabel('|\psi_L(t)|');
subplot(3,1,3); plot(f,abs(Psi_L)); xlabel('f'); ylabel('|\Psi_L(f)|');

% Ptychography performance
figure(4)
subplot(2,2,1); plot(Fid); title('Fid');
subplot(2,2,2); plot(log10(Ea)); title('log10 Ea');
subplot(2,2,3); plot(log10(Df)); title('log10 Df');


