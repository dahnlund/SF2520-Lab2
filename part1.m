%%% Part 1 of CE2. David Ahnlund and Emil Gestsson
clear, clc;
%Coefficients
a = 0.1; b = 0.4; Q0 = 7000; alpha0 = 50; Tout = 25; T0 = 100;

%% A

v = 1;
N_list = [10 20 40 80 160 320 640 1280];
saved_T = zeros(length(N_list),1);

for i = 1:length(N_list)

N = N_list(i);
h = 1/N;

z = 0:h:1;


D1 = (-1/h^2 - v/(2*h)) * ones(N-2,1);
D2 = (2/h^2) * ones(N-1,1);
D3 = (v/(2*h)-1/h^2) * ones(N-2,1);
A = diag(D1, -1) + diag(D2, 0) + diag(D3, 1);

Q_func = @(z) Q0 * sin((z-a)*pi / (b-a)) .* (a<=z).*(z<=b);
f = Q_func(z(2:end-1))';

% For the boundary T(0)
%--------------------------------|
f(1) = f(1) - (-1/h^2 - v/(2*h))*T0;
%--------------------------------|

% For the boundary T(1):
%--------------------------------|
alpha = @(v) sqrt(v^2/4 + alpha0^2) - v/2;

theta = (v/(2*h) - 1/h^2) / (3/(2*h)+alpha(v));

A(end,end-1) = (-1/h^2 - v/(2*h)-theta/(2*h));
A(end,end) = (2/h^2 + 4/(2*h)*theta);
f(end) = f(end) - theta * alpha(v)*Tout;
%--------------------------------|

T = A\f;
T = [T0; T];

saved_T(i) = T(z==0.9);
plot(z(1:end-1),T)
hold on

if ismember(N, [80 160 320])
    fprintf("T at z = 0.5 for N = %.0f: %.3f\n", N,T(z==0.5))
end
end

% Check convergence rate (should be 2)
d = diff(saved_T);
errors = sqrt(sum(d.^2, 2));
diff(log2(flip(errors)))

%% Part B


v_list = [1 5 15 100];


for i = 1:length(v_list)

N = 500;
h = 1/N;

z = 0:h:1;

v = v_list(i);

D1 = (-1/h^2 - v/(2*h)) * ones(N-2,1);
D2 = (2/h^2) * ones(N-1,1);
D3 = (v/(2*h)-1/h^2) * ones(N-2,1);
A = diag(D1, -1) + diag(D2, 0) + diag(D3, 1);

Q_func = @(z) Q0 * sin((z-a)*pi / (b-a)) .* (a<=z).*(z<=b);
f = Q_func(z(2:end-1))';

% For the boundary T(0)
%--------------------------------|
f(1) = f(1) - (-1/h^2 - v/(2*h))*T0;
%--------------------------------|

% For the boundary T(1):
%--------------------------------|
alpha = @(v) sqrt(v^2/4 + alpha0^2) - v/2;

theta = (v/(2*h) - 1/h^2) / (3/(2*h)+alpha(v));

A(end,end-1) = (-1/h^2 - v/(2*h)-theta/(2*h));
A(end,end) = (2/h^2 + 4/(2*h)*theta);
f(end) = f(end) - theta * alpha(v)*Tout;
%--------------------------------|

T = A\f;
T = [T0; T];

plot(z(1:end-1),T)
hold on


end
