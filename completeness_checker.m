%% Ptycho Completeness Checker
% Verifica numericamente se um conjunto de observáveis pticográficos é completo. Esses operadores
% serão da forma
%       O_kla = P_l F_a' |k><k| F_a P_l ,
% sendo que os P_l são projetores em subespaçcos do espaço de Hilbert de interesse, |k> são os
% estados da base computacional, F_a é uma transf. de Fourier fracional de ordem a e F_a' a sua
% adjunta.
%
% Esse script usa s função frft para as transformadas de Fourier fracionais.



%% Parametros
% Caso interessante pra mais tarde: usar mais de um valor em r, a e proj._family pra verificar a completeza (pensar
% nessas variáveis como vetores/cell ?)

d   = 10;                   % Dimensão do espaço de Hilbert
r   = 6;                    % Rank dos projetores, normalmente por volta de d/2
a   = [0.3]                 % ordens da frFT; a*pi/2 é o angulo da rotação da representação de Wigner;
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

% Transformadas fracionais de  Fourier
Na  = length(a);                    % núm. de ordens da frFT a serem usadas

foo = zeros(d,1);                   % var. auxiliar com os estados da base computacional (escritos na propria base comp.)
foo(1) = 1;

frF = zeros(d,d,Na);                 % frF(:,:,a) = matriz de F_a
for s = 1:d
    for t=1:Na
        frF(:,s) = frft(foo,a(t));       % aplica a frft no estado da base foo e salva o resultado na coluna de Fa
        foo = circshift(foo,[1,0]);     % shift pra montar o póximo estado da base
    end
end

% Aplicar nos estados da base comp. e salvar o resultado
foo = zeros(d,d);                   % var. auxiliar com os estados da base computacional (escritos na propria base comp.)
foo(1,1) = 1;

O = zeros(d,d,d,Nproj,Na);          % O(:,:,k,l,a) = O_kla
for s=1:d
    for t=1:Nproj
        for u=1:Na
            O(:,:,s,t,u)  = P(:,:,t)*frF(:,:,u)'*foo*frF(:,:,u)*P(:,:,t);
        end
    end
    foo = circshift(foo,[1,1]);
end




%% Verificando a completeza

if(Na~=1)
    warning('mais de um valor de a está sendo usado... o teste de rank pode não ser recomendado')
end

M   = zeros(d^2,d^2);
for s=1:Na
    for t=1:Nproj
        for u=1:d
            ind         = u + (t-1)*d + (s-1)*Nproj*d;
            M(:,ind)	= reshape(O(:,:,u,t,s).',[d^2 1]);
        end
    end
end

% Estratégia 1: montar a matriz expandida a partir dos observaveis O e calcular seu rank; isso só deve
% funcionar se estivermos usando um valor individual para a
RM	= rank(M);
disp(['Rank =	' num2str(RM)])
if(RM < d^2)
    disp('The set of observables is rank deficient.')
else
    disp('The set of observables is full rank.')
end



% Estratégia 2: encontrar autovalores e descartar os muito pequenos
lam	= eig(M);
thresh	= 0.1/d^2;
nlam	= sum( abs(lam)<thresh );

disp(['Number of small eigenvalues:	' num2str(nlam)])
disp(['Smallest eigenvalue modulum:	' num2str(abs(lam(end)))])



