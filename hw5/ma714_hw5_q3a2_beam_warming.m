h = 0.025; % gap of x
k = 0.015; % gap of t
T = 0.6;

% h = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
% k = 1e-3;
% T = 0.1;
%%%%%%%%%%%%%%%
n = T / k + 1; % number of t
m = 1 / h + 1; % number of x

x = linspace(0, 1, m);
t = linspace(0, T, n);
%%%%%%%%%%%%%%%
e = ones(m, 1);
A = spdiags([1 * e, -4 * e, 3 * e, 1 * e, -4 * e], [-2, -1, 0, m - 2, m - 1], m, m); % generate tridiagonal matrix A
B = spdiags([1 * e, -2 * e, 1 * e, 1 * e, -2 * e], [-2, -1, 0, m - 2, m - 1], m, m); % generate tridiagonal matrix B

U = initial_2(x);
U_store = U;

for i = 1:(n - 1)
    U = U - k / (2 * h) * A * U + k * k / (2 * h * h) * B * U;
    U_store = [U_store, U];
end

[xs, ts] = meshgrid(x, t);
U_store = reshape(U_store, [n, n]);

surf(xs, ts, U_store, 'facealpha', 0.5);
% % norm(U - ext, inf)

% xlim ([0, 1])
% ylim ([0, 1])
% zlim auto
% set(gca, 'XDir', 'reverse')
% set(gca, 'YDir', 'reverse')
xlabel("x");
ylabel("t");
zlabel("U");
title("The second initial condition, Beam-Warming, h = " + h + ", k = " + k);

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
