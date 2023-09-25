%%% Part 2 of CE2. David Ahnlund and Emil Gestsson
clear, clc;
%Coefficients
Lx = 12; Ly = 5; T_ext = 25;

%% A

N = 60; M = 25;
h = Lx/N;  %should be same as Ly/M
f = 2 * ones(N,M);

Sx = -1/h^2 * (diag(-ones(N-2,1),-1)+diag(2*ones(N-1,1),0) + diag(-ones(N-2,1),1));
Sy = -1/h^2 * (diag(-ones(M-2,1),-1)+diag(2*ones(M-1,1),0) + diag(-ones(M-2,1),1));

A = kron(eye(size(Sy)),Sx) + kron(Sy, eye(size(Sx)));



