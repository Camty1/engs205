%% ENGS 205 Chebyshev Derivative Test
% Cameron Wolfe

%% Setup
u_analytical = @(t, x) exp(-pi^2 * t) .* exp(sin(pi * x));
dudx0_analytical = @(x) pi * exp(sin(pi * x)) .* cos(pi * x);

x = -1:.02:1;

u0 = u_analytical(0, x);
dudx0 = dudx0_analytical(x);

%% Chebyshev calculation
figure(1);
plot(x, u0, 'DisplayName', 'Analytical');
hold on;

figure(2);
plot(x, dudx0, 'DisplayName', 'Analytical');
hold on;

N_vals = 3:64;
plot_N_vals = [4 8 16 32 64];
max_error = [];
rms_error = [];
deriv_max_error = [];
deriv_rms_error = [];
for N=N_vals
    % Calculate Chebyshev polynomials
    chebyshev_polys = {[1], [1 0]};
    for i=2:N
        chebyshev_polys{i+1} = 2*[chebyshev_polys{i} 0] - [0 0 chebyshev_polys{i-1}];
    end

    % Chebyshev coefficients
    u_gauss_lobatto = u_analytical(0, cos(pi * (0:N) / N));
    u_gauss_lobatto = [u_gauss_lobatto u_gauss_lobatto(end-1:-1:2)];

    u_hat = fft(u_gauss_lobatto);
    u_hat = real(u_hat(1:N+1));
    u_hat = u_hat / N ./ resize(c(0:N, N), size(u_hat));
    u_hat_1 = zeros(size(u_hat));
    u_hat_1(end-1) = 2 * N * u_hat(end);

    for k=N-2:-1:0
        u_hat_1(k+1) = u_hat_1(k+3) + 2 * (k+1) * u_hat(k+2);
    end

    u_hat_1(1) = u_hat_1(1) / 2;

    % Interpolation
    u_cheby = zeros(size(x));
    u_prime_cheby = zeros(size(x));

    for i=1:N+1
        u_cheby = u_cheby + real(u_hat(i)) * polyval(chebyshev_polys{i}, x);
        u_prime_cheby = u_prime_cheby + real(u_hat_1(i)) * polyval(chebyshev_polys{i}, x);
    end

    % Error calculation
    error = u_cheby - u0;
    max_error(end+1) = max(abs(error));
    rms_error(end+1) = sqrt(mean(error.^2));
    deriv_error = u_prime_cheby - dudx0;
    deriv_max_error(end+1) = max(abs(deriv_error));
    deriv_rms_error(end+1) = sqrt(mean(deriv_error.^2));

    if ismember(N, plot_N_vals)
        figure(1);
        plot(x, u_cheby, '--', 'DisplayName', "N = " + num2str(N));
        figure(2);
        plot(x, u_prime_cheby, '--', 'DisplayName', "N = " + num2str(N));
    end

end

figure(1);
hold off;
title("Chebyshev Approximations of u(x) by Order (N)");
xlabel("x");
ylabel("y");
legend();

figure(2);
hold off;
title("Chebyshev Approximations of u'(x) by Order (N)");
xlabel("x");
ylabel("y");
legend();

figure(3);
loglog(N_vals, max_error);
hold on;
loglog(N_vals, rms_error);
hold off;
title("Error of Chebyshev Approximations of u(x) vs Order");
xlabel("Order (N)");
ylabel("Error");
legend(["Max Error", "RMS Error"]);

figure(4);
loglog(N_vals, deriv_max_error);
hold on;
loglog(N_vals, deriv_rms_error);
hold off;
title("Error of Chebyshev Approximations of u'(x) vs Order");
xlabel("Order (N)");
ylabel("Error");
legend(["Max Error", "RMS Error"]);

function out = c(k, N) 
    out = ones(size(k));
    out(k == 0 | k == N) = 2;
end
