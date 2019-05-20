function [Obj_fin,Ea,Df,Fid,bad_res] =  spectralPIE_prop(par_ptycho,par_prop,f,Obj_ini,S,A_frF,Psi)

% Obj_ini should be the initial envelope+phase guess, in the freq. domain
% Psi will be the correct envelope+phase
% the frequency axis will be centered at zero,but not necessarily simetrical (fftshift thing...)

% figure(1); plot(1:length(Psi),ifftshift(A_frF))
% figure(3); plot(1:length(Psi),A_frF,1:length(Psi),abs(Psi),'o')

bad_res     = 0;
stp_flag    = 0;
npty    = 0;
Ea  = [];
Df  = [];
Fid = [];
Obj = Obj_ini;

switch par_ptycho.mom_mode 
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
    switch par_ptycho.cycling_method
        case 'linear'
            s = 1:par_ptycho.Nmasks;
        case 'random'
            s = randperm(par_ptycho.Nmasks);
        otherwise
            error('Unrecognized cycling method')
    end
    
    % Ptychographic iteration:
    % -(we begin in the Fourier plane)
    % -apply freq. filter
    % -propagate forth to intermediate plane
    % -correct ampitudes (A_frF)
    % -propagate to Fourier plane (apply forward frFT)
    % -correct estimate (Wiener) -> in the Fourier plane... that's where I know the illumination
    
    %f0  = C/par_wf.lamb0;
    %f0  = 0;
    td  = par_prop.L/par_prop.v;
    b   = par_prop.Dv*par_prop.L/pi;
    %H   = atten*exp(-1i*( 2*pi*td*(f-f0) + pi^2*b*(f-f0) ));
    H   = exp(-1i*( pi^2*b*f.^2 ));
    for r=1:par_ptycho.Nmasks
        G = S(s(r),:).*Obj;
        gfr = ifft( H.*G ); 
        gfr_c = A_frF(s(r),:).*exp(1i*angle(gfr));
        G_c = conj(H).*fft(gfr_c);
        
        U = updweight(S(s(r),:),par_ptycho.delta);
        Obj_new = Obj + par_ptycho.beta*U.*( G_c-G );
    end
    
    % Momentum update
    if(strcmp(par_ptycho.mom_mode,'classic')||strcmp(par_ptycho.mom_mode,'nesterov'))
        mom_t = mom_t+1;
        if(mom_t==par_ptycho.mom_T)
            v = par_ptycho.eta*v + Obj_new-Obj_old;
            switch par_ptycho.mom_mode
                case 'classic'
                    Obj_new =  Obj_old + v;
                case 'nesterov'
                    Obj_new = Obj_new + par_ptycho.eta*v;
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
    switch par_ptycho.stp_crit
        case 'rel-change'
            if(df<par_ptycho.dthresh || npty>=par_ptycho.npty_max)
                stp_flag = 1;
                if(npty>=par_ptycho.npty_max)
                    bad_res = 1;
                end
            end
        case 'npty'
            if(npty>=par_ptycho.npty_max)
                stp_flag = 1;
                if(df>par_ptycho.dthresh)
                    bad_res = 1;
                end
            end
        otherwise
            error('Unrecognized stop criterion')
    end
    
end

Obj_fin = Obj;