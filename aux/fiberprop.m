function [psiL,H] = fiberprop(lamb0,t,psi,L,Dv,v,atten)
    % psi   -> pulse wavefunction
    % lamb0 -> central wavelength of psi
    % t     -> time axis of psi
    % L     -> length of the fiber
    % Dv    -> Group Veocity Dispersion of the fiber
    % v     -> group velocity
    % atten -> attenuation of the fiber
    % Fundamentas of photonics, Saleh and Teich, ch. 22.3 (p. 962)
    % I will use the transfer function in order to simuate the propagation.
    % https://www.rp-photonics.com/group_velocity_dispersion.html -> GVD of fibers

    C = 3e+08;       %(speed of light in SI)

    N = length(t);
    dt = abs(t(2)-t(1));
    Fs = 1/dt;
    f = (0:N-1)*Fs/N;
    f = fftshift(f-mean(f));
    
%     f0 = C/lamb0;
    f0 = 0;

    td = L/v;
    b = Dv*L/pi;

%     H = atten*exp(-1i*( 2*pi*td*(f-f0) + pi^2*b*(f-f0) ));
    H = atten*exp(-1i*( pi^2*b*(f-f0).^2 ));

    Psi = fft(psi);
    PsiL = Psi.*H;
    psiL = ifft(PsiL);
end



