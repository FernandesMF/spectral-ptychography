

% cleaning yalmip memory
yalmip('clear');

F = class('double');

% defining the SDP variables
Rho = sdpvar(df, df, 'hermitian', 'complex');

% standard constraints
F = [rho_dash>=0,trace(Rho)==1];

% observables

Obs=cell(Nobs,df,df);


Prob = size(Nobs);

Delta = sdpvar(Nobs,'Real');

F = [F,Delta>=0];

for i=1:Nobs
    measure = trace(Rho*Obs{i});
    F = [F,measure>=Prob(i)(1-Delta(i))];
    F = [F,measure<=Prob(i)(1-Delta(i))];
end



% cost function
E = sdpvar(1,1,'Real');
F = [F,E>=0];
E = sum(Delta);



SOLUTION = solvesdp(F,E);
disp('DEBUGGING');
problema = doulbe(SOLUTION.problem);
disp(yalmiperror(problema));

Rho = double(Rho);
Delta = double(Delta);