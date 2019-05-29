clear all

% Setting parameters
d = 4;      % Hilbert space dimension
r = 2;      % rank of mask projectors; for now, we will use the contiguous family with d projectors scheme
a = [0.5 1];     % orders of frft that will be "measured" (1 --> )

Nproj = d;
Nobs = d*Nproj*numel(a);     % ptychographic projectors * comp. basis projectors

% Picking a random state
rho = cubitt_RandomDensityMatrix(d);

% Cleaning yalmip memory
yalmip('clear');
F = class('double');

% Defining the SDP variables
Rho = sdpvar(d,d,'hermitian','complex');
Delta = sdpvar(Nobs,1);

% Standard constraints
F = [Rho>=0,trace(Rho)==1];
F = [F,Delta>=0];

% Observables
Obs = zeros(d,d,Nobs);
pro = zeros(1,d);
pro(1:r) = 1;
com = zeros(1,d);
com(1) = 1;

w = exp(-1i*2*pi/d);
% jk = kron((0:d-1),(0:d-1)'); --> wrong matrix
jk = zeros(d,d);
for j=1:d-1
    for k=1:d-1
        jk(j+1,k+1) = j*k;
    end
end
Fou = w.^jk/sqrt(d);
iFou = (w).'.^jk/sqrt(d);

for s=1:d
    for l = 1:Nproj
        Pl = diag(circshift(pro,[0,l-1]));
        Com = circshift(com,[0,s-1]).'*circshift(com,[0,s-1]);
        foo = (s-1)*Nproj + l;
        Obs(:,:,foo) = Pl*iFou*Com*Fou*Pl;
    end
end

% Measured probabilities
Prob = zeros(Nobs,1);

for m = 1:Nobs
    Prob = trace(Rho*Obs(:,:,m));
    measured(m) = real(trace(rho*Obs(:,:,m)));
    F = [F,Prob-measured(m)*(1-Delta(m))>=0];
    F = [F,Prob-measured(m)*(1+Delta(m))<=0];
end

% Cost function
% here the cost function is the sum of the "flexibilizations" around the measured values
E = sdpvar(1,1);
F = [F,E>=0];
E = sum(Delta);

% Solving
SOLUTION = solvesdp(F,E,sdpsettings('solve','mosek'));
disp('DEBUGGING');
problema = double(SOLUTION.problem);
disp(yalmiperror(problema));

Rho = double(Rho)
Delta = double(Delta.')
dist = hsDistance(rho,Rho)
fid = Fidelity(rho,Rho)