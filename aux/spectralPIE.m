function [Obj_fin,Ea,Df,Fid,bad_res] =  spectralPIE(par_ptycho,Obj_ini,S,A_frF,Psi)

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
    % -apply inverse frFT (propagate back to intermediate plane)
    % -correct ampitudes (A_frF)
    % -propagate to Fourier plane (apply forward frFT)
    % -correct estimate (Wiener) -> in the Fourier plane... that's where I know the illumination
    for r=1:par_ptycho.Nmasks
        G = S(s(r),:).*Obj;
        gfr = frft(G,-1+par_ptycho.frft_a).';
        gfr_c = A_frF(s(r),:).*exp(1i*angle(gfr));
        G_c = frft(gfr_c,1-par_ptycho.frft_a).';    % i guess this is the troubled line...
        
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