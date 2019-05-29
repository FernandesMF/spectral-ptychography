%% Ptycho Completeness Checker
% Levanta numericamente o menor autovalor dos operadores de MUBs.

clear all


%% Parametros

d   = 5;                    % Dimensão do espaço de Hilbert
tol = 1e-14;                % Tolerância/epsilon em torno do zero (teste de ortogonalidade)
test    = 1;                % Fazer testes de orgotonalidade e viés?
ign	= 6;                    % Índice da base a ser ignorada


%% Montando a matriz com vetores das mubs

% Checando a dimensão
if(~isprime(d))
    error('Dimension d is not prime.')
end

% Gerando os estados das MUBs
%   Receita do artigo de Wooters e Fields: 
%   Optimal state-determination by mutually unbiased measurements, Annals of Physics 1991 (1989), p. 370

vecs	= zeros(d,d,d+1);

%   base comp.
vecs(:,:,1)	= diag(ones(d,1));

%   outras mubs
%       v_rkl = exp( 1i*2*pi*(rl²+kl)/d )/sqrt(d)    (eq 11)
%           r - índice da base; r = 1, ... , d
%           k - índice do vetor da base; k = 0, ... , d-1
%           l - componente do vetor (na base computacional); l = 0, ... , d-1
for r=1:d
    for k=0:d-1
        for l=0:d-1
            vecs(l+1,k+1,r+1)	= exp( 1i*2*pi*(r*l^2 + k*l)/d )/sqrt(d);
        end
    end
end

% Conferindo a falta de viés
if(test)
%   teste de ortogonalidade dentro de cada bases
    flag	= 0;
    for r=1:d+1
        for k=1:d-1
            for l=k+1:d
                foo	= vecs(:,k,r)'*vecs(:,l,r);
                if( foo > tol )
                    warning('Falha no teste de ortogonaildade:\t <%d,%d|%d,%d> = %f',k,r,l,r,foo);
                    flag    = 1;
                end
            end
        end
    end
    if(~flag)
        sprintf('Teste de ortogonalidade ok.')
    end

%   teste de falta de viés entre as bases
    flag    = 0;
    for r1=1:d
        for r2=r1+1:d+1
            for k=1:d
                for l=1:d
                    foo	= abs(vecs(:,k,r1)'*vecs(:,l,r2));
                    if( abs(foo-1/sqrt(d)) > tol )
                        warning('Falha no teste de viés:\t <%d,%d|%d,%d> = %f',k,r1,l,r2,foo);
                        flag    = 1;
                    end
                end
            end
        end        
    end
    if(~flag)
        sprintf('Teste de viés ok.')
    end
    
end     %if(test)



%% Montando matriz do esquema de medidas e calculando autovalores
%   Vou ignorar uma das bases na hora de construir (variável ign)

M   = zeros(d^2,d^2);
ind = 1:d+1; ind(ign) = [];
s   = 1;
for r=ind
    for k=1:d
        P	= vecs(:,k,r)*vecs(:,k,r)';
        M(:,s)	= reshape(P.',[d^2 1]);
        s   = s+1;
    end
end
lam	= abs(eig(M))
lammin	= min(lam)










