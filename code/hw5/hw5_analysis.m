u_true = readmatrix("problem_definition/uf_1024.dat");
x_true = linspace(0, 2*pi, 1025);

N_vals = [8, 16, 32, 64, 128, 256, 512, 1024];

max_error = zeros(2, length(N_vals));
rms_error = zeros(2, length(N_vals));

close all;
for i=1:length(N_vals)
    N = N_vals(i);
    ua = readmatrix(['problem_definition/uf_' num2str(N) '.dat']);
    up = readmatrix(['output/uf_pseudo_' num2str(N) '.dat']);
    uc = readmatrix(['output/uf_coll_' num2str(N) '.dat']);

    error = up - ua;
    max_error(1, i) = max(abs(error));
    rms_error(1, i) = sqrt(mean(error.^2));

    error = uc - ua;
    max_error(2, i) = max(abs(error));
    rms_error(2, i) = sqrt(mean(error.^2));

    figure();
    plot(x_true, u_true([1:end 1]), 'HandleVisibility', 'off');
    hold on;
    plot(2 * pi / N * (0:(N-1)), up, '+', 'DisplayName', 'Pseudospectral');
    plot(2 * pi / N * (0:(N-1)), uc, 'x', 'DisplayName', 'Collocation');
    hold off;
    title("Burger's Equation Spectral Solution N = " + num2str(N)); 
    xlabel("x");
    ylabel("y");
    xlim([0, 2 * pi]);
    legend(location="northwest");
end

figure();
loglog(N_vals, max_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, max_error(2, :), '--x', 'DisplayName', 'Collocation');
hold off;
title("Burger's Equations Spectral Solution Max Error vs N");
xlabel("N");
ylabel("Max Error");
legend(location="northeast");


figure();
loglog(N_vals, rms_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, rms_error(2, :), '--x', 'DisplayName', 'Collocation');
hold off;
title("Burger's Equations Spectral Solution RMS Error vs N");
xlabel("N");
ylabel("RMS Error");
legend(location="northeast");

u_true = readmatrix("problem_definition/uf_pi_4_1024.dat");
max_error = zeros(2, length(N_vals));
rms_error = zeros(2, length(N_vals));

for i=1:length(N_vals)
    N = N_vals(i);
    ua = readmatrix(['problem_definition/uf_pi_4_' num2str(N) '.dat']);
    up = readmatrix(['output/uf_pseudo_pi_4_' num2str(N) '.dat']);
    uc = readmatrix(['output/uf_coll_pi_4_' num2str(N) '.dat']);

    error = up - ua;
    max_error(1, i) = max(abs(error));
    rms_error(1, i) = sqrt(mean(error.^2));

    error = uc - ua;
    max_error(2, i) = max(abs(error));
    rms_error(2, i) = sqrt(mean(error.^2));

    figure();
    plot(x_true, u_true([1:end 1]), 'HandleVisibility', 'off');
    hold on;
    plot(2 * pi / N * (0:(N-1)), up, '+', 'DisplayName', 'Pseudospectral');
    plot(2 * pi / N * (0:(N-1)), uc, 'x', 'DisplayName', 'Collocation');
    hold off;
    title("Burger's Equation Spectral Solution N = " + num2str(N)); 
    xlabel("x");
    ylabel("y");
    xlim([0, 2 * pi]);
    legend(location="northwest");
end

figure();
loglog(N_vals, max_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, max_error(2, :), '--x', 'DisplayName', 'Collocation');
hold off;
title("Burger's Equations Spectral Solution Max Error vs N");
xlabel("N");
ylabel("Max Error");
legend(location="northeast");

figure();
loglog(N_vals, rms_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, rms_error(2, :), '--x', 'DisplayName', 'Collocation');
hold off;
title("Burger's Equations Spectral Solution RMS Error vs N");
xlabel("N");
ylabel("RMS Error");
legend(location="northeast");


u_true = readmatrix("problem_definition/uf_1024.dat");
anti_max_error = zeros(1, length(N_vals))
anti_rms_error = zeros(1, length(N_vals))

for i=1:length(N_vals)
    N = N_vals(i);
    ut = readmatrix(['problem_definition/uf_' num2str(N) '.dat']);
    ua = readmatrix(['output/uf_anti_' num2str(N) '.dat']);

    error = ua - ut;
    anti_max_error(1, i) = max(abs(error));
    anti_rms_error(1, i) = sqrt(mean(error.^2));

    figure();
    plot(x_true, u_true([1:end 1]), 'HandleVisibility', 'off');
    hold on;
    plot(2 * pi / N * (0:(N-1)), ua, '+', 'DisplayName', 'Antialiased');
    hold off;
    title("Burger's Equation Spectral Solution N = " + num2str(N)); 
    xlabel("x");
    ylabel("y");
    xlim([0, 2 * pi]);
    legend(location="northwest");
end

figure();
loglog(N_vals, max_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, max_error(2, :), '--x', 'DisplayName', 'Collocation');
loglog(N_vals, anti_max_error, '-.o', 'DisplayName', 'Antialiased');
hold off;
title("Burger's Equations Spectral Solution Max Error vs N");
xlabel("N");
ylabel("Max Error");
legend(location="northeast");

figure();
loglog(N_vals, rms_error(1,:), '-+', 'DisplayName', 'Pseudospectral');
hold on;
loglog(N_vals, rms_error(2, :), '--x', 'DisplayName', 'Collocation');
loglog(N_vals, anti_rms_error, '-.o', 'DisplayName', 'Antialiased');
hold off;
title("Burger's Equations Spectral Solution RMS Error vs N");
xlabel("N");
ylabel("RMS Error");
legend(location="northeast");
