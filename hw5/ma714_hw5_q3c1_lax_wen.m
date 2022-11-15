% H = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
H = [1e-2, 2e-3, 4e-4];
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

    U = initial_1(x);
    U_true = true_1(x, T);
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
% U
% plot(x, U);
% hold on;
% plot(x, U_true);
% hold off;
% xlabel("x");
% ylabel("U");
% legend("U", "U_true")

% function output = initial_1(x)
%     output = 2 .* exp(- (100 .* x - 50).^2/10);
%     output = output';
% end

% function output = initial_2(x)
%     output = exp(-100 .* (x - 0.5).^2) .* sin(80 .* pi .* x);
%     output = output';
% end

% function output = true_1(x, t)
%     output = 2 .* exp(- (100 .* (mod((x - t), 1)) - 50).^2/10);
%     output = output';
% end

% function output = true_2(x)
%     output = exp(-100 .* ((mod((x - t), 1)) - 0.5).^2) .* sin(80 .* pi .* (mod((x - t), 1)));
%     output = output';
% end

% end

% hold off;
% ma714_hw5_q3b_beam_warming(0.01); hold on;
% ma714_hw5_q3b_beam_warming(0.02);
% ma714_hw5_q3b_beam_warming(0.05);
% ma714_hw5_q3b_beam_warming(0.1);
% legend("0.01", "0.02", "0.05", "0.1");
% title("The second initial condition, Beam-Warming, T = 1, h = 0.01, 0.02, 0.05, 0.1, k = 0.01");

% err1 = ma714_hw5_q3b_beam_warming(0.01);
% err2 = ma714_hw5_q3b_beam_warming(0.02);
% err3 = ma714_hw5_q3b_beam_warming(0.05);
% err4 = ma714_hw5_q3b_beam_warming(0.1);
% err = [err1, err2, err3, err4];
% h = [0.01, 0.02, 0.05, 0.1];
% [slope, intercept] = er_order(h, err);
title("Log-Log Error with the first initial condition using Lax Wendroff, slope = " + slope);
