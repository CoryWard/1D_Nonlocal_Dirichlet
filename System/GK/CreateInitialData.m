%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the par structure containing the parameters and 
% system definition. This file is essentially a contructor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%%%%%
% Create Computational Domain
%%%%%

par.Lw = 20; % The computational domain is [-Lw, Lw]
par.M = 200; % Number of equally spaced points between [0,Lw] not including 0.
             % Also note that M has to be EVEN!
par.h = par.Lw/par.M; % Space between points
par.L = par.Lw/2; % BVP is defined on [-L,L]

par.x = par.h*(-par.M:1:par.M)'; % Computational Domain
par.xb = par.h*(-par.M/2+1:par.M/2-1)'; % BVP Domain
par.xbe = par.h*(-par.M/2:par.M/2)'; % BVP Domain, including the boundary points

%%%%%
% Define Kernel
%%%%%

nu = @(y) (1/2)*exp(-abs(y));   % Kernel
Fi = @(y) (1/2)*exp(-abs(y));    % Define second antiderivative
Fip = @(y) -(1/2)*sign(y).*exp(-abs(y)); % Define first antiderivative

par.F = Fi(par.x);
par.Fp = Fip(par.x);

%%%%%
% Define Forcing Function
%%%%%

fi = @(y) sech(y) - (1/2)*exp(-y).*log(1+exp(2*y))-...
        (1/2)*exp(y).*log(1+exp(2*y)) + y.*exp(y);

par.f = fi(par.xb);

%%%%%
% Calculate Corresponding Integrals
%%%%%

par.f1 = (1/2)*(1/par.h^2)*(2 - exp(-par.h)*(par.h^2 + 2*par.h + 2));
par.A = exp(-par.Lw);

%%%%%
% Define Dirichlet BC and corresponding Integral
%%%%%

g = @(y) sech(y);
Bi = @(x) (1/2)*exp(x).*log(exp(-2*par.Lw)+exp(2*x)) - exp(x).*x ...
        + (1/2)*exp(-x).*log(exp(-2*par.Lw)+exp(-2*x)) + exp(-x).*x;
par.B = Bi(par.xb);

%%%%%
% Define the omega_i
%%%%%

omega(1) = par.f1 - Fip(par.h) +(Fi(2*par.h) -Fi(par.h))/par.h;
omega(2:par.M-1,1) = (1/par.h)*(Fi(par.h*(2:par.M-1)' + par.h)-...
                    2*Fi(par.h*(2:par.M-1)') + Fi(par.h*(2:par.M-1)' - par.h));
omega(par.M) = Fip(par.h*par.M) + (Fi((par.M-1)*par.h) - Fi(par.M*par.h))/par.h;

par.omega = [omega(par.M:-1:1); 0; omega];

%%%%%
% Define the c_i
%%%%%

par.c1 = sum(par.omega);

for i = -par.M/2+1:par.M/2-1          % There's a better way to do this.
    par.c2(i+par.M/2,1)=0;                
    for j = -par.M:i-par.M/2
        par.c2(i+par.M/2) = par.c2(i+par.M/2) + g((i-j)*par.h).*par.omega(j+par.M+1);
    end
end

for i = -par.M/2+1:par.M/2 -1         
    par.c3(i+par.M/2,1)=0;                
    for j = i+par.M/2:par.M
        par.c3(i+par.M/2) = par.c3(i+par.M/2) + g((i-j)*par.h).*par.omega(j+par.M+1);
    end
end

%%%%%
% Define the matrix w_hat
%%%%%

par.wh = toeplitz(par.omega(par.M+1:end-2),par.omega(par.M+1:-1:3));

%%%%%
% Create path to save data
%%%%%

filename='GK1';
file_names.Data = strcat('System/GK/Data/',filename,'.mat');

save(strcat('Data/',filename,'.mat'), 'par', 'file_names');
