%%% Part 2 of CE2. David Ahnlund and Emil Gestsson
clear, clc;
%Coefficients
Lx = 12; Ly = 5; T_ext = 25;

%% A

N = 60*2; M = 25*2;
h = Lx/N;  %should be same as Ly/M
F = 2 * ones(N-1,M-1);

Sx = -1/h^2 * (diag(-ones(N-2,1),-1)+diag(2*ones(N-1,1),0) + diag(-ones(N-2,1),1));
Sy = -1/h^2 * (diag(-ones(M-2,1),-1)+diag(2*ones(M-1,1),0) + diag(-ones(M-2,1),1));

%Boundary condition for x
Sx(1,1) = 2/(3*h^2); Sx(1,2) = -2/(3*h^2);
Sx(end,end) = 2/(3*h^2); Sx(end, end-1) = 4/(3*h^2);

%Boundary condition for y
Sy(end,end) = 2/(3*h^2); Sy(end, end-1) = 4/(3*h^2);

A = kron(eye(size(Sy)),Sx) + kron(Sy, eye(size(Sx)));

F(:,1) = F(:,1) - T_ext;

f = reshape(F, (N-1)*(M-1),1);

u = A\f;

U = reshape(u, (N-1), (M-1));

mesh(U)