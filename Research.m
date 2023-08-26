% This is a MATLAB code for 1D Finite Difference Method
% This code explains the solutions with uniform mesh and Shishkin meshes
% The model problems is singularly perturbed convection diffusion problem in 1D
% with Dirithlet boundary conditions
% -epsilon*u" + u' = f with u(0)=0, u(1)=0
% This code also explains the maximum norm error

close
clc
format short
epsilon=10^-6;  % Singularly perturbed parameter
p=[32 64 128 256 512 1024 2048 4096 8192]; % Generate different meshes
for j=1:length(p) % a for loop to run through all the meshes are record errors and rates

%% Generate a uniform mesh
n=p(j); % Number of interior points
h=1/(n+1); % Step size
y=0:h:1; % The interval (the domain)
x=y';

%% Generate a Shishkin mesh
% beta=0.9;
% sigma=2;
% N=p(j); % number of nodel point in the domain (an even number)
% tau = min(1/2, sigma*epsilon*log(N)/beta);
%     x = unique([linspace(0,1-tau, N/2), ...
%         linspace(1-tau,1, N/2)])';
% n = length(x)-2;
% h=zeros(n,1);
% for i = 1:n-1
% h(i) = x(i+1) - x(i); % step size in the Shishkin mesh
% end
% h(end)=h(end-1);
%% Build matrices to get Finite Difference solution
A =(2*diag(ones(1,n)) - 1*diag(ones(1,n-1),1) - 1*diag(ones(1,n-1),-1)); % Stiffness matrix
C =-0.5*h.*(diag(zeros(1,n)) - 1*diag(ones(1,n-1),1) + 1*diag(ones(1,n-1),-1));%Convection matrix
f=h.^2.*x(2:end-1); % Source term
u_int=(epsilon*A+C)\f; % Finite difference approximation for interior points
U=zeros(n+2,1); % Sparse matrix to save all the solutions
U(1)=0; % Left boundary condition
U(end)=0; % Right  boundary condition
U= [U(1); u_int; U(end)] ; % Complete numerical approximations
%% Exact solution
v1=x.*(0.5*x + epsilon);
v2=(0.5+ epsilon);
v3=exp((x-1)/epsilon) - exp(-1/epsilon);
v4=1 - exp(-1/epsilon);
u_exact= v1-v2.*(v3/v4);
%% record Maximum norm errors
evec=U-u_exact;
maxerr=max(evec)
if j>2
rate= abs(log(errold/maxerr)/log(2));
end
errold=maxerr;
end
%% Plot Solutions
figure(1)
plot(x,U,'r*',x,u_exact,'b','LineWidth',2);
legend('Numerical Solution', 'Exact Solution')
title('Solution Comparison')
