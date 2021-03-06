clear all
% close all

% ATTENTION: use (nano) SI units throughout this code
C = 3e+08;  % nm/ns, speed of light

%% Parameters

% Time and frequency vectors 
dt  = 1e-06;        % ns; always leave this as 1e-XX, or there will be problem when resampling.
T   = 1e-01;        % ns
t   = -T/2:dt:T/2;	% ns
N   = length(t);
Fs  = 1/dt;         % GHz
df  = Fs/N;         % GHz
f   = (0:N-1)*df;	% GHz
f   = f-mean(f);

% Photon gaussian wavefunction
% gaussian with quadratic phase profile (which can be pretty much any... this is controllable)
par_wf.a   = 0.1;               % chirp coeficient (no units)
par_wf.lam0    = 1560e-00;      % nm
par_wf.sig_wl  = 3e-00;         % nm

par_wf.f0  = C/par_wf.lam0;     % GHz
par_wf.sig_f   = par_wf.sig_wl*C/par_wf.lam0^2;     % GHz
par_wf.sig_t   = 1/(2*pi*par_wf.sig_f);             % ns

% Interference filter
par_if.R    = 0.99;          % reflectance
par_if.T    = 0.01;          % transmisttance
par_if.dr   = 0;             % phase gain at reflection (radians)
par_if.dt   = 0;             % phase gain at transmission (radians)
par_if.n    = 1.45;          % refraction index
par_if.d    = 6.5e+03;       % nm; distance between interfaces
par_if.thetas  = [0:1:10]*pi/180.0;             % in radians
par_if.F    = 4*par_if.R/( 1-par_if.R )^2;      % finesse

% Propagation through fiber
par_prop.L  = 100;     % m
par_prop.Dl = 18e-06;   % s/m^2
par_prop.Dv = -par_wf.lam0^2/(2*pi*C)*par_prop.Dl;
par_prop.v  = 1.5e+08;


% Detection
% par_detec.timeres       = 1.115e-03;	%ns; estimated in the script detec_timeres; make equal to dt if resampling unwanted
par_detec.timeres       = 1.115e-03/5;	%ns; estimated in the script detec_timeres; make equal to dt if resampling unwanted
par_detec.flag_noise	= 0;            % 0: no noise; 1: only jitter; 2: only poissonian; 3: jitter and poissonian
par_detec.av_counts     = 1e+04;
par_detec.jitter        = 5e-03;       %ns



% Ptychography algorithm parameters
par_ptycho.init_est_method  = 'random';
par_ptycho.cycling_method   = 'random';  %'random', 'linear' 
par_ptycho.stp_crit     = 'npty';        %'rel-change', 'npty'
par_ptycho.dthresh      = 1e-02;         % either the rel change threshold or the number of ptycho. iterations
par_ptycho.npty_max     = 1000;
par_ptycho.n_rerand     = 100;          
par_ptycho.beta     = 0.5;               % step size/correction multiplier
par_ptycho.delta    = 0.05;              % Wiener filter/rPIE parameter
par_ptycho.Nmasks   = length(par_if.thetas);    % number of spectral masks to be used

par_ptycho.mom_mode     = 'nesterov';    %'none','classic','nesterov'
par_ptycho.mom_T    = 30;                % number of iterations to gather when using momentum
par_ptycho.eta      = 0.9;               % momentum 'friction' (or better yet, '1-friction')

par_ptycho.frft_a   = 2*atan(4*pi^2*par_prop.Dv*par_prop.L/(par_detec.timeres^2))/pi; %FIX: change dt to detection dt




%% Generate photon initial state
psi = exp(-t.^2/(2*par_wf.sig_t^2) ).*exp(-1i*2*pi*par_wf.f0.*t).*...
       exp(1i*par_wf.a*t.^2/par_wf.sig_t^2 );
% psi = exp(-t.^2/(2*par_wf.sig_t^2) ).*exp(1i*par_wf.a*t.^2/par_wf.sig_t^2 );
psi = fftshift(psi);
t   = fftshift(t);
f   = fftshift(f);
Psi = fft(psi);


%% Generate spectral masks
I = zeros(length(par_if.thetas),N);
I = ones(size(I)); % for debugging purposes...
for s=1:length(par_if.thetas)
    delta 	= par_if.dr + pi*f*par_if.n*par_if.d*cos( par_if.thetas(s) )/C;
%     I(s,:)  = par_if.T^2/( (1-par_if.R)^2 ) * 1./( 1+par_if.F*(sin( delta )).^2 );    %intensity mask...
%     I(s,:)  = par_if.T*exp(1i*2*par_if.dt)./( 1 - par_if.R*exp(1i*2*par_if.dr)*exp(1i*delta) );
end



%% Propagate filtered wavefunctions
Psi_filt    = zeros(length(par_if.thetas),N);
Psi_L       = zeros(length(par_if.thetas),N);
psi_L       = zeros(length(par_if.thetas),N);
for s=1:length(par_if.thetas)
    Psi_filt(s,:)   = I(s,:).*Psi;
    %[Psi_L(s,:),H]  = fiberprop(par_wf.lam0,t,Psi_filt(s,:),par_prop.L,par_prop.Dv,par_prop.v,1);
    
    %f0  = C/par_wf.lamb0;
    f0  = 0;
    td  = par_prop.L/par_prop.v;
    b   = par_prop.Dv*par_prop.L/pi;
    %H   = atten*exp(-1i*( 2*pi*td*(f-f0) + pi^2*b*(f-f0) ));
    H   = exp(-1i*( pi^2*b*(f-f0).^2 ));
    Psi_L(s,:)  = Psi_filt(s,:).*H;
    psi_L(s,:)  = ifft(Psi_L(s,:));
end




%% Generate measurements

% Applying jitter
if(par_detec.flag_noise==1||par_detec.flag_noise==3)
    disp('apply jitter')
end

% Resampling signals at detector resolution
if(par_detec.timeres<dt)
    warning('Detector time resolution is set to smaller than original dt... are you sure about this?')
end

[rsp,rsq]   = rat(dt/par_detec.timeres,1e-08);
if(rsp~=1)
    error('Resampling wont work, rsp is not 1; did you set dt different than 1e-XX?')
end

% rs_psi = psi(1:rsq:end);
rs_t   = t(1:rsq:end);
rs_N   = length(rs_t);
rs_dt  = abs(rs_t(2)-rs_t(1));

par_ptycho.frft_a   = 2*atan(4*pi^2*par_prop.Dv*par_prop.L/(rs_dt^2))/pi %FIX: change dt to detection dt

rs_Fs  = 1/rs_dt;           % GHz
rs_df  = rs_Fs/rs_N;        % GHz
rs_f   = (0:rs_N-1)*rs_df;	% GHz
rs_f   = rs_f-mean(rs_f)-par_wf.f0; % GHz

if( ~issorted(rs_f) || ~issorted(ifftshift(f)) )
    error('Unsorted frequency vector; this will spoil the resampling')
end

foo	= zeros(size(rs_f));
fooo = zeros(size(rs_t));
for s=1:length(foo)
    [~,foo(s)]	= min( (f-rs_f(s)).^2 );
    [~,fooo(s)] = min( (t-rs_t(s)).^2 );
end
% rs_psi = psi(fooo);
% rs_Psi = fft(rs_psi);
% rs_Psi = Psi(fooo);

rs_psi = exp(-rs_t.^2/(2*par_wf.sig_t^2) ).*exp(-1i*2*pi*par_wf.f0.*rs_t).*...
       exp(1i*par_wf.a*rs_t.^2/par_wf.sig_t^2 );
% psi = exp(-t.^2/(2*par_wf.sig_t^2) ).*exp(1i*par_wf.a*t.^2/par_wf.sig_t^2 );
% rs_psi = fftshift(rs_psi);
% t   = fftshift(t);
% f   = fftshift(f);
rs_Psi = fft(ifftshift(rs_psi));

% FIX: get wf/I at the points newly defined in the scales;
for s=1:length(par_if.thetas)
    rs_psi_L(s,:)	= psi_L(s,1:rsq:end);
    rs_I(s,:)	= I(s,foo);
end

% Checking resampling
figure(21)
subplot(2,2,1); plot(1:rs_N,rs_t,'-o',(1:N)*rs_N/N,t); title('rs\_t');
subplot(2,2,2); plot(1:rs_N,rs_f,'-o',(1:N)*rs_N/N,f); title('rs\_f');
subplot(2,2,3); plot(rs_t,abs(rs_psi),'-o',t,abs(psi)); title('\psi(t)'); legend('resampled','original')
subplot(2,2,4); plot(rs_f,abs(rs_Psi),'-o',f,abs(Psi)); title('\Psi(f)'); legend('resampled','original')

figure(22)
subplot(2,1,1); plot(rs_t,abs(rs_psi_L)); title('resampled \psi_{L}(t)');
subplot(2,1,2); plot(t,abs(psi_L)); title('original \psi_{L}(t)');

figure(23)
subplot(3,1,1); plot(rs_f,abs(rs_I)); title('resampled I(f)');
subplot(3,1,2); plot(rs_f,abs(Psi(foo))); title('resampled \Psi(f)');
subplot(3,1,3); plot(f,abs(I)); title('original I(f)');

% Applying poissonian noise
if(par_detec.flag_noise==2||par_detec.flag_noise==3)
    disp('apply poissonian noise')
end



%% Run ptychography
%FIX: implement rerandomization
disp('Starting Ptychography')

% Initial random estimate (smoothed, gaussian complex numbers)
% obj     = (randn(1,rs_N)-0.5)+1i*(randn(1,rs_N)-0.5);
obj = rs_psi;           % sanity test
obj_ini	= obj;          % saving initial estimate, just for checking purposes...
% obj = smooth(obj,10).';


A	= sqrt(sum( abs(psi).^2 ));
A_obj   = sqrt(sum( abs(obj).^2 ));
obj = A*obj/A_obj;      % normalizing initial estimate
Obj = fft(obj);
Obj = rs_Psi;           % sanity test
Obj_ini = Obj;          % saving initial estimate for checking purposes

% Ptychographic engine
[Obj,Ea,Df,Fid,bad_res] =  spectralPIE(par_ptycho,Obj,rs_I,abs(rs_psi_L),rs_Psi);





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

% Original and filtered wavefunctions (time and freq. domains)
figure(3)
subplot(3,1,1); plot(t,abs(psi)); xlabel('t'); ylabel('|\psi(t)|');
subplot(3,1,2); plot(t,abs(psi_L)); xlabel('t'); ylabel('|\psi_L(t)|');
subplot(3,1,3); plot(f,abs(Psi_L)); xlabel('f'); ylabel('|\Psi_L(f)|');

% Ptychography performance
figure(4)
subplot(2,2,1); plot(Fid); title('Fid');
subplot(2,2,2); plot(log10(Ea)); title('log10 Ea');
subplot(2,2,3); plot(log10(Df)); title('log10 Df');

% Sanity test
figure(5)
plot(rs_f,abs(Obj_ini),f,abs(Psi)); legend('Obj_ini','\Psi(f)')



%% Old code

% psi = exp(-t.^2/(2*sig_t^2) ).*exp(1i*2*pi*f0.*t).*exp(1i*a*t.^2/sig_t^2 );
% psi = fftshift(psi);
% t   = fftshift(t);
% [psiL] = fiberprop(lam0,t,psi,L,Dv,v,atten)



% from https://en.wikipedia.org/wiki/Wave_packet , section about dispersive solution; I will use x=0
% further

% psi = 1./sqrt(1+1i*2*t).*exp( -(k0*t).^2./(1+4*t.^2) ).*exp( 1i*(k0^2*t)./(4+16*t.^2) );
% psi = ifftshift(psi);
% Psi = fft(psi);

%{
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

figure(1)
subplot(1,2,1)
plot(t,abs(psi));
subplot(1,2,2)
plot(t,angle(psi))

%}