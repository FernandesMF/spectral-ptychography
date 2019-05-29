function a = annihilfock_mobral(N)

% Essa rotina faz um operador destruição super simples.

% N é o número de Fock máximo. Exemplo, se N = 4, esse
% aniquilador mata até 4 fótons e tá no espaço de até
% quatro fótons


n = 1:N;
a = diag(sqrt(n),1); 