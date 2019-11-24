clear all
% close all

% FIX: implement rerandomization

% ATTENTION: use (nano) SI units throughout this code
C = 3e+08;  % nm/ns, speed of light

%% Parameters

% Execution parameters
chck_resample   = 0;        % Make plots to check resampling?
chck_jitter     = 0;        % Plots to check jitter blurring?
chck_arr_times  = 1;        % Plots to check generated arrival times?

% Time and frequency vectors 
dt  = 1e-06;        % ns; always leave this as 1e-XX, or there will be problem when resampling.
T   = 5e-02;        % ns
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
par_prop.L  = 1;     % m
par_prop.Dl = 18e-06;   % s/m^2
par_prop.Dv = -par_wf.lam0^2/(2*pi*C)*par_prop.Dl;
par_prop.v  = 1.5e+08;


% Detection
% par_detec.timeres       = 1.115e-03;	%ns; estimated in the script detec_timeres; make equal to dt if resampling unwanted
par_detec.timeres       = 1.115e-03/5;	%ns; estimated in the script detec_timeres; make equal to dt if resampling unwanted
par_detec.flag_noise	= 2;            % 0: no noise; 1: only jitter; 2: only poissonian; 3: jitter and poissonian
par_detec.av_counts     = 1e+06;        
par_detec.jitter        = 5e-03;        %ns; good jitter, from possible detector in the future
% par_detec.jitter        = 25e-03;       %ns; bad jitter, from current detectors



% Ptychography algorithm parameters
par_ptycho.init_est_method  = 'random';
par_ptycho.cycling_method   = 'linear';  %'random', 'linear' 
par_ptycho.stp_crit     = 'npty';        %'rel-change', 'npty'
par_ptycho.dthresh      = 1e-02;         % either the rel change threshold or the number of ptycho. iterations
par_ptycho.npty_max     = 100;
par_ptycho.n_rerand     = 100;          
par_ptycho.beta     = 0.5;               % step size/correction multiplier
par_ptycho.delta    = 0.05;              % Wiener filter/rPIE parameter
par_ptycho.Nmasks   = length(par_if.thetas);    % number of spectral masks to be used

par_ptycho.mom_mode     = 'none';    %'none','classic','nesterov'
par_ptycho.mom_T    = 30;                % number of iterations to gather when using momentum
par_ptycho.eta      = 0.9;               % momentum 'friction' (or better yet, '1-friction')

par_ptycho.frft_a   = [];
% par_ptycho.frft_a   = 2*atan(4*pi^2*par_prop.Dv*par_prop.L/(par_detec.timeres^2))/pi; %FIX: change dt to detection dt




%% Generate photon initial state
% psi = exp(-t.^2/(2*par_wf.sig_t^2) ).*exp(-1i*2*pi*par_wf.f0.*t).*...
%        exp(1i*par_wf.a*t.^2/par_wf.sig_t^2 );
% psi = exp(-t.^2/(2*par_wf.sig_t^2) ).*exp(1i*par_wf.a*t.^2/par_wf.sig_t^2 );
psi = exp(-t.^2/(2*2*par_wf.sig_t^2) ).*exp(1i*par_wf.a*t.^2/(2*par_wf.sig_t^2) ); % sigs. of intensity, not wf!
psi = psi/sqrt(sum(abs(psi).^2));       %/(2*pi*2*par_wf.sig_t^2) did not work... normalizing wf (remember its std is sqrt(2)*sit]g_t)
psi = ifftshift(psi);
t   = ifftshift(t);
f   = fftshift(f);
Psi = fft(psi);
% WARNING: if the form of psi changes, I should correct 'ptycho_resample' function.


%% Generate spectral masks
I = zeros(length(par_if.thetas),N);
% I = ones(size(I)); % for debugging purposes...
for s=1:length(par_if.thetas)
    delta 	= par_if.dr + pi*(f+par_wf.f0)*par_if.n*par_if.d*cos( par_if.thetas(s) )/C;
%     I(s,:)  = par_if.T^2/( (1-par_if.R)^2 ) * 1./( 1+par_if.F*(sin( delta )).^2 );    %intensity mask...
    I(s,:)  = par_if.T*exp(1i*2*par_if.dt)./( 1 - par_if.R*exp(1i*2*par_if.dr)*exp(1i*delta) );
end
% FIX: check if I is increasing the total intensity of psi_L...



%% Propagate filtered wavefunctions
Psi_filt    = zeros(length(par_if.thetas),N);
Psi_L       = zeros(length(par_if.thetas),N);
psi_L       = zeros(length(par_if.thetas),N);
for s=1:length(par_if.thetas)
    Psi_filt(s,:)   = I(s,:).*Psi;
    %[Psi_L(s,:),H]  = fiberprop(par_wf.lam0,t,Psi_filt(s,:),par_prop.L,par_prop.Dv,par_prop.v,1);
    
    %f0  = C/par_wf.lamb0;
    %f0  = 0;
    td  = par_prop.L/par_prop.v;
    b   = par_prop.Dv*par_prop.L/pi;
    %H   = atten*exp(-1i*( 2*pi*td*(f-f0) + pi^2*b*(f-f0) ));
    H   = exp(-1i*( pi^2*b*f.^2 ));
    Psi_L(s,:)  = Psi_filt(s,:).*H;
    psi_L(s,:)  = ifft(Psi_L(s,:));
end




%% Generate measurements
% noise flag values:
%	-0: no noise
%   -1: only jitter
%   -2: only poissonian
%   -3: jitter and poissonian

% Calling resampling function
disp('Resampling');
[rs_t,rs_f,rs_N,rs_dt,rs_Fs,rs_df,rs_psi,rs_psi_L,rs_Psi,rs_I,par_ptycho] = ...
               ptycho_resample(par_ptycho,par_detec,par_wf,par_prop,par_if,t,dt,f,psi_L,I);
% Checking resampling
if(chck_resample)
    figure(21)
    subplot(2,2,1); plot(1:rs_N,rs_t,'-o',(1:N)*rs_N/N,t); title('rs\_t');
    subplot(2,2,2); plot(1:rs_N,rs_f,'-o',(1:N)*rs_N/N,f*rs_N/N); title('rs\_f (original f rescaled)');
    subplot(2,2,3); plot(rs_t,abs(rs_psi),'-o',t,abs(psi)); title('\psi(t)'); legend('resampled','original')
    subplot(2,2,4); plot(rs_f,abs(rs_Psi)*N/rs_N,'-o',f,abs(Psi)); title('\Psi(f)'); legend('resampled (rescaled)','original')

    figure(22)
    subplot(2,1,1); plot(rs_t,abs(rs_psi_L)); title('resampled \psi_{L}(t)');
    subplot(2,1,2); plot(t,abs(psi_L)); title('original \psi_{L}(t)');

    figure(23)
    subplot(3,1,1); plot(rs_f,abs(rs_I)); title('resampled I(f)');
    subplot(3,1,2); plot(rs_f,abs(rs_Psi)); title('resampled \Psi(f)');
    subplot(3,1,3); plot(f,abs(I)); title('original I(f)');
end

% Applying jitter to psi_L's (which are the detected quantities)
if(par_detec.flag_noise==1||par_detec.flag_noise==3)
    disp('Applying jitter');
    blur	= exp( -(t-mean(t)).^2/(2*par_detec.jitter^2) )/sqrt(2*pi*par_detec.jitter^2);
    foo     = psi_L;
    for s=1:length(par_if.thetas)
        psi_L(s,:)	= sqrt( conv( abs(psi_L(s,:)).^2,blur,'same' ) );
        psi_L(s,:)  = ifftshift(psi_L(s,:));
    end
    % Checking jitter
    if(chck_jitter)
        figure(31)
        subplot(3,1,1); plot(t,abs(foo)); title('original psi_L');
        subplot(3,1,2); plot(t,blur); title('blur');
        subplot(3,1,3); plot(t,abs(psi_L)); title('blurred psi_L');
    end
end

% Generating counts
 if(par_detec.flag_noise==2||par_detec.flag_noise==3)
    disp('Generating counts')
    % Generating Cumulative Distribution Functions (of each |psi_L|²)
    totI_L	= zeros(length(par_if.thetas),1);      % Total intensity of psi_L
    cdf_L	= zeros(size(psi_L));                   % CDF of |psi_L|²
    for s=1:length(par_if.thetas)
        foo = 0;        % accumulator variable to walk through psi_L
        for r=1:N
            foo	= foo + abs(psi_L(s,r)).^2;
            cdf_L(s,r)	= foo;
        end
        totI_L(s)	= foo;
    end
    
    % Generating number of counts (poissonian) for each CDF (prop. to total intensity)
    totCounts_L = poissrnd(par_detec.av_counts*totI_L);
%     totCounts_L = poissrnd(totI_L);
    % FIX: totI_L has many more counts than I expected... perhaps I_L is increasing amplitudes?
    
    % Drawing random arrival times from CDF
    sim_psi_L	= zeros(size(psi_L));   
    for s=1:length(par_if.thetas)
        foo	= totI_L(s)*rand(totCounts_L(s),1);
        foo = sort(foo,'ascend');
        bar = 1;
        dum = 0;
        for r=1:N
            while( bar<totCounts_L(s) && cdf_L(s,r)>foo(bar) )
                bar = bar+1;
            end
            sim_psi_L(s,r)	= bar-dum;
            dum = bar;
        end
    end
    
    % Check generated arrival times
    if(chck_arr_times)
        figure(32)
        subplot(2,1,1); plot(t,abs(psi_L).^2); title('original pdfs');
        subplot(2,1,2); plot(t,sim_psi_L); title('simulated data');
        figure(33)
        subplot(2,1,1); plot(t,abs(psi_L)); title('sqrt of original pdfs');
        subplot(2,1,2); plot(t,sqrt(sim_psi_L)); title('sqrt of simulated data');
    end
end



%% Run ptychography
%FIX: implement rerandomization
disp('Starting Ptychography')

% Initial random estimate (smoothed, gaussian complex numbers)
obj	= (randn(1,N)-0.5)+1i*(randn(1,N)-0.5);
obj	= smooth(obj,9).';
% obj = psi;            % sanity test
obj_ini	= obj;          % saving initial estimate, just for checking purposes...


A	= sqrt(sum( abs(psi).^2 ));
A_obj   = sqrt(sum( abs(obj).^2 ));
obj = A*obj/A_obj;      % normalizing initial estimate
Obj	= fft(obj);
% Obj = Psi;            % sanity test
Obj_ini	= Obj;          % saving initial estimate for checking purposes

% Ptychographic engine
% FIX: make a clause to choose which of these lines to call
% [Obj,Ea,Df,Fid,bad_res] =  spectralPIE(par_ptycho,Obj,rs_I,abs(rs_psi_L),rs_Psi);
% [Obj_fin,Ea,Df,Fid,bad_res] =  spectralPIE_prop(par_ptycho,par_prop,f,Obj_ini,S,A_frF,Psi)

% [Obj,Ea,Df,Fid,bad_res] =  spectralPIE_prop(par_ptycho,par_prop,par_detec,...
%                                             f,t,Obj,I,abs(psi_L),Psi,);
[Obj,Ea,Df,Fid,bad_res] =  spectralPIE_prop(par_ptycho,par_prop,par_detec,...
                                            f,t,Obj,I,sqrt(sim_psi_L),Psi);





%% Show results
disp('Plotting results');

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

% Original and recovered functions
figure(5)
plot(f,abs(Obj),f,abs(Psi)); legend('Obj','\Psi(f)')


