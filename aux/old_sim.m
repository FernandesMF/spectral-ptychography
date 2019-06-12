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

