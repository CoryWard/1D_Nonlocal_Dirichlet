function u = NonlocalSystemSolver(par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the nonlocal linear system 
%           Lu = b
% defined by the par structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% Set parameter values for ease of use.
%%%%%

c1 = par.c1;
c2 = par.c2;
c3 = par.c3;
wh = par.wh;
A = par.A;
B = par.B;
f = par.f;
M = par.M;

%%%%%
% Define system
%%%%%

L = c1*eye(M-1,M-1) - wh + A*eye(M-1,M-1);  % Matrix
b = c2 + c3 + f + B;

%%%%%
% Solve system
%%%%%

u = L\b;
end

