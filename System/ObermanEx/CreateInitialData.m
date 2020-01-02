%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the par structure containing the parameters and 
% system definition. This file is essentially a contructor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%%%%%
% Create Computational Domain
%%%%%

par.Lw = 2; % The computational domain is [-Lw, Lw]
par.M = 40; % Number of equally spaced points between [0,Lw] not including 0.
             % Also note that M has to be even
par.h = par.Lw/par.M; % Space between points
par.L = par.Lw/2; % BVP is defined on [-L,L]

par.x = par.h*(-par.M:1:par.M)'; % Computational Domain
par.xb = par.h*(-par.M/2+1:par.M/2-1)'; % BVP Domain
par.xbe = par.h*(-par.M/2:par.M/2)'; % BVP Domain, including the boundary points

%%%%%
% Define Kernel
%%%%%

par.a = .8; % a can't be one!
Ca = (par.a*2^(par.a-1)*gamma((par.a+1)/2))/(pi^(1/2)*gamma((2-par.a)/2));
nu = @(y) Ca*abs(y).^(-1-a);   % Kernel
Fi = @(y) Ca/((par.a-1)*par.a)*abs(y).^(1-par.a);    % Define second antiderivative
Fip = @(y) -Ca./(par.a).*sign(y).*abs(y).^(-par.a);   % Define first antiderivative

par.F = Fi(par.x);
par.Fp = Fip(par.x);

%%%%%
% Define Forcing Function
%%%%%

fi = @(y) 0.*y+1;
par.f = fi(par.xb);

%%%%%
% Calculate Corresponding Integrals
%%%%%

par.f1 = Ca/(2-par.a)*par.h^(-par.a);
par.A = 2*Ca/(par.a*par.Lw^par.a);

%%%%%
% Define Dirichlet BC and correspond Integral
%%%%%

g = @(y) 0.*y;
Bi = @(x) 0.*x;
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

for i = -par.M/2+1:par.M/2-1          
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

filename='ObEx1';
file_names.Data = strcat('System/ObermanEx/Data/',filename,'.mat');

save(strcat('Data/',filename,'.mat'), 'par','file_names');
