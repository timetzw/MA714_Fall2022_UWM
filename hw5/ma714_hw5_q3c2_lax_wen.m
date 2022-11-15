% H = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
z = [500, 1500, 4500];
% H = [2e-3, 1e-3, ];
H = 1 ./ z;
K = H / 20;
% k = 1e-3;
T = 0.05;

%%%%%%%%%%%%%%%
errs = []

for j = 1:length(H)
    h = H(j);
    k = K(j);
    n = T / k + 1; % number of t
    m = 1 / h + 1; % number of x

    x = linspace(0, 1, m);
    t = linspace(0, T, n);
    %%%%%%%%%%%%%%%
    e = ones(m, 1);
    A = spdiags([2 * e, 1 * e, -6 * e, 3 * e, 2 * e, 1 * e, -6 * e], [1 - m, -2, -1, 0, 1, m - 2, m - 1], m, m); % generate tridiagonal matrix A
    B = spdiags([1 * e, 1 * e, -2 * e, 1 * e, 1 * e], [1 - m, -1, 0, 1, m - 1], m, m); % generate tridiagonal matrix B
    C = spdiags([1 * e, -1 * e, 3 * e, -3 * e, 1 * e, -1 * e, 3 * e], [1 - m, -2, -1, 0, 1, m - 2, m - 1], m, m); % generate tridiagonal matrix C

    U = initial_2(x);
    U_true = true_2(x, T);
    % U = initial_2(x);
    % U_true = true_2(x,T);

    for i = 1:(n - 1)
        U = U - k / (6 * h) * A * U + k * k / (2 * h * h) * B * U - k * k * k / (6 * h * h * h) * C * U;
    end

    err = U - U_true;
    err = norm(err, inf);
    errs = [errs, err];
end

errs
[slope, intercept] = er_order(H, errs);
slope
title("Log-Log Error with the second initial condition using Lax Wendroff, slope = " + slope);
