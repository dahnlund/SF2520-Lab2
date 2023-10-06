%%% Part 1 of CE2. David Ahnlund and Emil Gestsson
clear, clc;
%Coefficients
a = 0.1; b = 0.4; Q0 = 7000; alpha0 = 50; Tout = 25; T0 = 100;

%% A

v = 1;
N_list = [10 20 40 80 160 320 640 1280 1280*2 1280*4];

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
T_N = 1/(3/(2*h)+alpha(v)) * (2*T(end)/h-T(end-1)/(2*h)+alpha(v)*Tout);
T = [T0; T; T_N];

saved_T(i) = T(z==1);
if ismember(N, [10 20 40 80])
    plot(z,T)
    hold on
end

if ismember(N, [80 160 320])
    fprintf("T(z = 0.5) = %.3f, for N = %.0f\n", T(z == 0.5), N)
end
end

for i = 2:length(N_list)-1
    fprintf("&(%.0f,%.0f)", N_list(i), 2*N_list(i))
end
fprintf("\n")

xlabel("Position throughout cylinder, z")
ylabel("Temperature, T")
legend("N = 10","N = 20","N = 40","N = 80")

% Check convergence rate (should be 2)
d = abs(diff(saved_T));
approxconv = abs(diff(log2(d)));

%Create table
fprintf("\n\n\n")
disp("Convergence rate analysis:")
for i = 2:length(N_list)-1
    fprintf("$(%.0f,%.0f)$&",N_list(i),N_list(i+1))
end
fprintf("\n")

for i = 1:length(approxconv)
    fprintf("$%.05f$&",approxconv(i))

end
fprintf("\n\n\n")

%% Part B


v_list = [1 5 15 100];
N = 1000;
h = 1/N;

for i = 1:length(v_list)

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
T_N = 1/(3/(2*h)+alpha(v)) * (2*T(end)/h-T(end-1)/(2*h)+alpha(v)*Tout);
T = [T0; T; T_N];

fprintf("T(z = 0.5) = %.3f, for v = %.0f\n", T(z == 0.5), v)

plot(z,T)
hold on

xlabel("Position throughout cylinder, z")
ylabel("Temperature, T")
legend("v = 1","v = 5","v = 15","v = 100")

end
