function [rs_t,rs_f,rs_N,rs_dt,rs_Fs,rs_df,rs_psi,rs_psi_L,rs_Psi,rs_I,par_ptycho] = ...
               ptycho_resample(par_ptycho,par_detec,par_wf,par_prop,par_if,t,dt,f,psi_L,I)
    % This function was previously part of the main script, but I made a function of it since it got
    % a bit extense.
    %
    % This function resamples:
    %   -the time and frequency scales
    %   -the time and frequency wavefunctions
    %   -the propagated time wavefunctions
    %   -the frequency masks.
    %
    % If psi is not a gaussian anymore, I should change the code here.
    
    

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

    par_ptycho.frft_a   = 2*atan(4*pi^2*par_prop.Dv*par_prop.L/(rs_dt^2))/pi; %FIX: change dt to detection dt

    rs_Fs  = 1/rs_dt;           % GHz
    rs_df  = rs_Fs/rs_N;        % GHz
    rs_f   = (0:rs_N-1)*rs_df;	% GHz
    rs_f   = rs_f-mean(rs_f);   % GHz

    if( ~issorted(rs_f) || ~issorted(ifftshift(f)) )
        error('Unsorted frequency vector; this will spoil the resampling')
    end

    foo	= zeros(size(rs_f));
    fooo	= zeros(size(rs_t));
    for s=1:length(foo)
        [~,foo(s)]	= min( (f-rs_f(s)).^2 );
        [~,fooo(s)] = min( (t-rs_t(s)).^2 );
    end
    % rs_psi = psi(fooo);
    % rs_Psi = fft(rs_psi);
    % rs_Psi = Psi(fooo);

    rs_psi	= exp(-rs_t.^2/(2*par_wf.sig_t^2) ).*exp(1i*par_wf.a*rs_t.^2/par_wf.sig_t^2 );
    % rs_psi = ifftshift(rs_psi);
    % rs_t   = ifftshift(rs_t);
    rs_f   = fftshift(rs_f);
    rs_Psi = fft(rs_psi);

    % FIX: get wf/I at the points newly defined in the scales;
    for s=1:length(par_if.thetas)
        rs_psi_L(s,:)	= psi_L(s,1:rsq:end);
    %     rs_psi_L(s,:)   = ifftshift(rs_psi_L(s,:));
        rs_I(s,:)	= I(s,foo);
        rs_I(s,:)   = fftshift(rs_I(s,:));
    end



end