%% Ptycho Completeness Checker
% Levanta numericamente o menor autovalor do conjunto de observadores pticográficos, que são da
% forma
%       O_kla = P_l F_a' |k><k| F_a P_l ,
% sendo que os P_l são projetores em subespaçcos do espaço de Hilbert de interesse, |k> são os
% estados da base computacional, F_a é uma transf. de Fourier fracional de ordem a e F_a' a sua
% adjunta.
%
% Esse script usa s função frft para as transformadas de Fourier fracionais.

clear all


%% Parametros
% Caso interessante pra mais tarde: usar mais de um valor em r, a e proj._family pra verificar a completeza (pensar
% nessas variáveis como vetores/cell ?)

d   = 10;                   % Dimensão do espaço de Hilbert
r   = 7;                    % Rank dos projetores, normalmente por volta de d/2
Na  = 1e03;                 % Número de valores de a
a   = linspace(0,1,Na);     % ordens da frFT; a*pi/2 é o angulo da rotação da representação de Wigner;
                            % nesse script, dá pra usar um valor individual ou vários valores pra a



%% Montando as matrizes dos operadores

% Projetores (P_l)
Nproj       = d;                    % núm. de projetores (=d agora, mas depois...)

pro         = zeros(1,d);           % variável auxiliar com a diagonal do primeiro Pl
pro(1:r)    = 1;

P           = zeros(d,d,Nproj);     % P(:,:,l) = matriz de Pl 
for s=1:Nproj
    P(:,:,s)    = diag(pro);                % coloca pro na diagonal de P
    pro         = circshift(pro,[0,1]);     % shift pra montar a diagonal do próximo projetor
end

lammin	= zeros(Na,1);
for z=1:Na
    % Transformadas fracionais de  Fourier
    foo = zeros(d,1);                   % var. auxiliar com os estados da base computacional (escritos na propria base comp.)
    foo(1) = 1;

    frFa = zeros(d,d);                 % frF(:,:) = matriz de F_a
    for s = 1:d
        frFa(:,s) = frft(foo,a(z));       % aplica a frft no estado da base foo e salva o resultado na coluna de Fa
        foo = circshift(foo,[1,0]);      % shift pra montar o póximo estado da base
    end

    % Aplicar nos estados da base comp. e salvar o resultado
    foo = zeros(d,d);                   % var. auxiliar com os estados da base computacional (escritos na propria base comp.)
    foo(1,1) = 1;

    O = zeros(d,d,d,Nproj);          % O(:,:,k,l,a) = O_kla
    for s=1:d
        for t=1:Nproj
            O(:,:,s,t)  = P(:,:,t)*frFa(:,:)'*foo*frFa(:,:)*P(:,:,t);
        end
        foo = circshift(foo,[1,1]);
    end
    
    % Montando matriz com operadores "vetorizados"
    M   = zeros(d^2,d^2);
    for t=1:Nproj
        for u=1:d
            ind         = u + (t-1)*d;
            M(:,ind)	= reshape(O(:,:,u,t).',[d^2 1]);
        end
    end
    lam	= abs(eig(M));
    lammin(z)	= min(lam);
end



%% Plots

fig = figure(1);
h	= plot(a,lammin);
xlabel('a')
ylabel('|\lambda_{\rm min}|')
title(['d = ' num2str(d) ', r = ' num2str(r)])

Width	= 800;
Height  = 600;
foo     = fig.Position;
fig. Position   = [foo(1) foo(2) Width Height];
h.LineWidth     = 1.5
h.Parent.FontSize   = 16

print(['./minabseig/minabseig_d' num2str(d) '_r' num2str(r) '.png'],'-dpng')





