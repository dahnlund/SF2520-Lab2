%%% Part 1 of CE2. David Ahnlund and Emil Gestsson

%Coefficients
a = 0.1; b = 0.4; Q0 = 7000; alpha0 = 50; Tout = 25; T0 = 100;

%% A

v = 1;
N = 100;
h = 1/N;

z = 0:h:1;


D1 = (-1/h^2 - v/(2*h)) * ones(N,1);
D2 = (2/h^2) * ones(N+1,1);
D3 = (v/(2*h)-1/h^2) * ones(N,1);
A = diag(D1, -1) + diag(D2, 0) + diag(D3, 1);

Q_func = @(z) Q0 * sin((z-a)*pi / (b-a)) .* (a<=z).*(z<=b);
Q = Q_func(z)';
Q(1) = Q(1) - (-1/h^2 - v/(2*h))*T0;  %Boundary condition at T(0)
T = A\Q;