function [X, J] = HowtoSDP(A)

%   author: Jessica Bavaresco
%   date: September, 2015                           

%   Minimizes Tr(XA) such that X is a quantum state.

%   INPUTS: 2x2 matrix A
%   OUTPUTS: optimal state X and the minimum of the function J = (Tr(XA))


% cleaning yalmip memory
yalmip('clear');

% defining the SDP variables
X = sdpvar(2, 2, 'hermitian', 'complex');

% standard constraints
F = [X>=0,trace(X)==1];

% cost function
J = trace(X*A);

% minimize J by varying X subject to F
SOLUTION = solvesdp(F,J,sdpsettings('solve','mosek','verbose',0))

%???verbose???,0: hides details of computation

X=double(X);
J=double(J);
