function U = updweight(P,delta)
    
    Pmax = max(abs(P));
%     U = abs(P).*conj(P)./((abs(P).^2+delta)*Pmax);    %Wiener filter
%     U = conj(P)./(Pmax^2+delta);                      %PIE
    U = conj(P)./( (1-delta)*abs(P).^2 + delta*Pmax^2 );%rPIE
    
end