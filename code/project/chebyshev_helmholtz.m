x = -1:.02:1;
u = 4 / (3 * pi^2) * (sin(pi * x / 2) - cos(pi * x)) + sin(pi * x / 2);
u_1 = 4 / pi^2 + (4 + pi^2) / pi^2 * sin(pi * x / 2);
u_simple = sin(pi * x / 2);

f = @(x) cos(pi * x);

figure(1);
plot(x, u, 'DisplayName', 'Analytical');
hold on;

N_vals = 3:64;
plot_N_vals = [4 8 16 32 64];
max_error = [];
rms_error = [];
deriv_max_error = [];
deriv_rms_error = [];
for N=N_vals
    chebyshev_polys = {[1], [1 0]};

    for i=2:N
        chebyshev_polys{i+1} = 2*[chebyshev_polys{i} 0] - [0 0 chebyshev_polys{i-1}];
    end

    D = zeros(N+1);
    for i=1:2:N
        D = D + 2 * diag(i:N, i);
    end

    A = D^2 + (pi / 2)^2 * eye(N+1);

    f_gauss_lobatto = f(cos(pi * (0:N) / N));
    f_gauss_lobatto = [f_gauss_lobatto f_gauss_lobatto(end-1:-1:2)];

    f_hat = fft(f_gauss_lobatto);
    f_hat = real(f_hat(1:N+1));
    f_hat = (f_hat / N ./ resize(c(0:N, N), size(f_hat)))';
    disp(f_hat)

    A(end-1, :) = (-1).^(0:N);
    A(end, :) = (-1).^(0:N) .* (0:N).^2;

    f_hat(end-1) = -1;
    f_hat(end) = 0;

    u_hat = A \ f_hat;

    u_cheby = zeros(size(x));

    for i=1:N+1
        u_cheby = u_cheby + real(u_hat(i)) * polyval(chebyshev_polys{i}, x);
    end

    error = u_cheby - u;
    max_error(end+1) = max(abs(error));
    rms_error(end+1) = sqrt(mean(error.^2));

    if ismember(N, plot_N_vals)
        figure(1);
        plot(x, u_cheby, '--', 'DisplayName', "N = " + num2str(N));
    end

end

figure(1);
hold off;
title("Chebyshev Approximations of u(x) by Order (N)");
xlabel("x");
ylabel("y");
legend(location="northwest");

figure(3);
loglog(N_vals, max_error);
hold on;
loglog(N_vals, rms_error);
hold off;
title("Error of Chebyshev Approximations of u(x) vs Order");
xlabel("Order (N)");
ylabel("Error");
legend(["Max Error", "RMS Error"]);

function out = c(k, N) 
    out = ones(size(k));
    out(k == 0 | k == N) = 2;
end