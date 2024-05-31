N = 8;
K = 4.0;
c = 4.0;
v = 0.2;
tf = pi/8;

dt = 8 / v / N^2;
m = pi / 8.0 / dt;
timesteps = ceil(m) * 100;
dt = pi / 8 / timesteps;

u = readmatrix(['problem_definition/u0_' num2str(N) '.dat']);
u_true = readmatrix(['problem_definition/uf_' num2str(N) '.dat']);
u_tilde = fft(u);

for l=1:timesteps
    u_prime_tilde = zeros(size(u_tilde));
    u_dbl_prime_tilde = zeros(size(u_tilde));
    for k=1:N
        if k <= N/2
            u_prime_tilde(k) = 1i * (k - 1) * u_tilde(k);
            u_dbl_prime_tilde(k) = (1i * (k - 1))^2 * u_tilde(k);
        else
            u_prime_tilde(k) = 1i * (k - N - 1) * u_tilde(k);
            u_dbl_prime_tilde(k) = (1i * (k - N - 1))^2 * u_tilde(k);
        end
    end

    u = ifft(u_tilde);
    u_prime = ifft(u_prime_tilde);

    u_u_prime = u .* u_prime;

    u_u_prime_tilde = fft(u_u_prime);

    u_tilde = u_tilde + dt * (v * u_dbl_prime_tilde - u_u_prime_tilde);
end

u = ifft(u_tilde);

plot(u_true);
hold on;
plot(real(u))
hold off;
