%%% Part 2 of CE2. David Ahnlund and Emil Gestsson
clear, clc;
%Coefficients
Lx = 12; Ly = 5; T_ext = 25;

%% A

N = 60; M = 25;
h = Lx/N;  %should be same as Ly/M
F = 2 * ones(N-1,M-1);

Sx = 1/h^2 * (diag(-ones(N-2,1),-1)+diag(2*ones(N-1,1),0) + diag(-ones(N-2,1),1));
Sy = 1/h^2 * (diag(-ones(M-2,1),-1)+diag(2*ones(M-1,1),0) + diag(-ones(M-2,1),1));

%Boundary condition for x
Sx(1,1) = 2/(3*h^2); Sx(1,2) = -2/(3*h^2);
Sx(end,end) = 2/(3*h^2); Sx(end, end-1) = -2/(3*h^2);

%Boundary condition for y
Sy(end,end) = 2/(3*h^2); Sy(end, end-1) = -2/(3*h^2);

A = kron(eye(size(Sy)),Sx) + kron(Sy, eye(size(Sx)));
A = sparse(A);

F(:,1) = F(:,1) + T_ext/h^2;

f = reshape(F, (N-1)*(M-1),1);

t = A\f;

T = reshape(t, (N-1), (M-1));
T0 = T_ext * ones(N-1,1);
T = [T0 T];

mesh(T)
xlabel("y")
ylabel("x")

x = 0:h:Lx-h;
y = 0:h:Ly-h;

fprintf("T(6,2) = %.3f, for N = %.0f\n", T(round(x,2)==6,round(y,2)==2), N)

%% b

N = 60;
h = Lx/N;
M = Ly/h;
f = 2;

%From anaylitcal solution we have:
%------------
c0 = T_ext;
c1 = f * Ly;
c2 = -f/2;
%------------


x = 0:h:Lx-h;
y = 0:h:Ly-h;

T_analytical = @(x,y) (c0 + c1*y + c2*y.^2).*ones(length(x),1);

mesh(T_analytical(x',y))
hold on
mesh(T)
xlabel("y")
ylabel("x")

fprintf("T(6,2) = %.3f, for analytical solution\n", T_analytical(6,2))


%% d

N = 60*4;
h = Lx/N;
M = Ly/h;

F_func = @(x,y) 100*exp(-1/2 * (x-4).^2 - 4*(y-1).^2);

x = h:h:Lx-h;
y = h:h:Ly-h;

F = F_func(x',y);

Sx = 1/h^2 * (diag(-ones(N-2,1),-1)+diag(2*ones(N-1,1),0) + diag(-ones(N-2,1),1));
Sy = 1/h^2 * (diag(-ones(M-2,1),-1)+diag(2*ones(M-1,1),0) + diag(-ones(M-2,1),1));

%Boundary condition for x
Sx(1,1) = 2/(3*h^2); Sx(1,2) = -2/(3*h^2);
Sx(end,end) = 2/(3*h^2); Sx(end, end-1) = -2/(3*h^2);

%Boundary condition for y
Sy(end,end) = 2/(3*h^2); Sy(end, end-1) = -2/(3*h^2);

A = kron(eye(size(Sy)),Sx) + kron(Sy, eye(size(Sx)));

A = sparse(A);

F(:,1) = F(:,1) + T_ext/h^2;

f = reshape(F, (N-1)*(M-1),1);

t = A\f;

T = reshape(t, (N-1), (M-1));
T0 = T_ext * ones(N-1,1);
T = [T0 T];

figure
mesh(T)
xlabel("y")
ylabel("x")

figure
imagesc(x,y,T)
figure
contour(T)

fprintf("T(6,2) = %.3f, for N = %.0f\n", T(round(x,6)==6,round(y,6)==2), N)
