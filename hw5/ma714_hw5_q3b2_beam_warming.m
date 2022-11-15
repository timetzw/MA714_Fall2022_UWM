H = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
k = 1e-3;
T = 0.1;
%%%%%%%%%%%%%%%

for j = 1:length(H)
    h = H(j);
    n = T / k + 1; % number of t
    m = 1 / h + 1; % number of x

    x = linspace(0, 1, m);
    t = linspace(0, T, n);
    %%%%%%%%%%%%%%%
    e = ones(m, 1);
    A = spdiags([1 * e, -4 * e, 3 * e, 1 * e, -4 * e], [-2, -1, 0, m - 2, m - 1], m, m); % generate tridiagonal matrix A
    B = spdiags([1 * e, -2 * e, 1 * e, 1 * e, -2 * e], [-2, -1, 0, m - 2, m - 1], m, m); % generate tridiagonal matrix B

    U = initial_2(x);

    for i = 1:(n - 1)
        U = U - k / (2 * h) * A * U + k * k / (2 * h * h) * B * U;
    end

    plot(x, U);
    hold on;
end

% U
% plot(x, U);
% hold on;
% plot(x, U_true);
% hold off;
xlabel("x");
ylabel("U");
legend("h=1e-1", "h=5e-2", "h=2.5e-2", "h=1.25e-2");
title("The second initial condition using Beam Warming at time T="+T);
