%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main:
% Solve the nonlocal linear system using the collocation
% method of Huang and Oberman. The file "CreateInitialData.m"
% acts as a constructor, defining the system, it's
% parameters, and the corresponding BVP. This file solves 
% the corresponding BVP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%
% Select which system to use
%%%%%

sys = load('System/GK/Data/GK1.mat');
% sys = load('System/ObermanEx/Data/ObEx1.mat');

par = sys.par;
file_names = sys.file_names;

%%%%%
% Solve system
%%%%%

u = NonlocalSystemSolver(par);

%%%%%
% Plot the solution.
%%%%%

%%% GK %%%
g=sech(par.xbe);
u = [sech(par.xbe(1)); u; sech(par.xbe(end))]; % Add the values of u at the 
                                                 % boundary points.

%%% Huang/Oberman %%%
% u=[0;u;0];
% g = (2^(-par.a)*gamma(1/2))/(gamma(1 +par.a/2)*...
%         gamma((1+par.a)/2))*(1-par.xbe.^2).^(par.a/2);


max(abs(u - g))
plot(par.xbe,u,'k',par.xbe,g,'r')

%%%%%
% Save solution to par structure
%%%%%

par.u = u;

save(file_names.Data,'par','file_names');