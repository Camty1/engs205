Ns = [4, 8, 16, 32, 64, 128, 256, 512, 1024];
K = 4;
c = 4;
v = 0.2;

burger_bounds = 50;

for N=Ns
    phi_eqn = @(x, t, n) exp(-(x - (2*n+1)*pi).^2 ./ (4 * v * t));
    dphi_eqn = @(x, t, n) -2 / (4 * v * t) .* exp(-(x - (2*n+1)*pi).^2 ./ (4 * v * t)) .* (x - (2 * n + 1)*pi);

    x = linspace(0, 2*pi, N+1);
    x = x(1:end-1)';
    phi = zeros(N, 1);
    dphi = zeros(N, 1);
    for n=-burger_bounds:burger_bounds
        phi = phi + phi_eqn(x + pi, 1, n);
        dphi = dphi + dphi_eqn(x + pi, 1, n);
    end
    u = K - 2 * v * dphi ./ phi;
    writematrix(u, ['problem_definition/u0_' num2str(N) '.dat']);

    phi = zeros(N, 1);
    dphi = zeros(N, 1);
    for n=-burger_bounds:burger_bounds
        phi = phi + phi_eqn(x - c * pi/8 + pi, pi/8 + 1, n);
        dphi = dphi + dphi_eqn(x - c * pi/8 + pi, pi/8 + 1, n);
    end
    u = K - 2 * v * dphi ./ phi;
    writematrix(u, ['problem_definition/uf_' num2str(N) '.dat']);

    phi = zeros(N, 1);
    dphi = zeros(N, 1);
    for n=-burger_bounds:burger_bounds
        phi = phi + phi_eqn(x - c * pi/4 + pi, pi/4 + 1, n);
        dphi = dphi + dphi_eqn(x - c * pi/4 + pi, pi/4 + 1, n);
    end
    u = K - 2 * v * dphi ./ phi;
    writematrix(u, "problem_definition/uf_pi_4_" + num2str(N) + ".dat");
end
