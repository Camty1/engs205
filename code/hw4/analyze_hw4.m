sizes = [4, 8, 16, 32, 64, 128, 256, 512, 1024];

close all;

dx = zeros(1, length(sizes));
error = zeros(4, length(sizes));
rmserror = zeros(size(error));

for i=1:length(sizes)
    x = linspace(0, 2*pi, sizes(i)+1);
    x = x(1:end-1);
    dx(i) = x(2) - x(1);

    analytical = readmatrix(['output/analytical_' num2str(sizes(i)) '.dat']);
    analytical_plot = readmatrix(['output/analytical_1024.dat']);
    fd = readmatrix(['output/fd_' num2str(sizes(i)) '.dat']);
    fem = readmatrix(['output/fem_' num2str(sizes(i)) '.dat']);
    sm = readmatrix(['output/sm_' num2str(sizes(i)) '.dat']);

    error(1, i) = max(abs(fd - analytical));
    error(2, i) = max(abs(fem - analytical));
    error(3, i) = max(abs(sm - analytical));

    rmserror(1, i) = sqrt(mean((fd - analytical).^2));
    rmserror(2, i) = sqrt(mean((fem - analytical).^2));
    rmserror(3, i) = sqrt(mean((sm - analytical).^2));

    figure(i);
    plot(linspace(0, 2*pi, 1024), analytical_plot, 'DisplayName', 'Analytical');
    hold on;
    plot([x, 2 * pi], fd([1:end 1]), 'DisplayName', 'Finite Difference');
    plot([x, 2 * pi], fem([1:end 1]), 'DisplayName', 'Finite Element Method');
    plot([x, 2 * pi], sm([1:end 1]), 'DisplayName', 'Spectral Method');
    hold off;
    title(['Approximation Comparisons for N = ' num2str(sizes(i))]);
    xlabel('x');
    ylabel("f'");
    xlim([0, 2*pi]);
    legend;
end

figure();
loglog(dx, error(1, :), 'DisplayName', 'Finite Difference');
hold on;
loglog(dx, error(2, :), 'DisplayName', 'Finite Element Method');
loglog(dx, error(3, :), 'DisplayName', 'Spectral Method');
hold off
legend('location', 'southwest');
title("Max Error of Approximation Methods vs Step Size");
xlabel("Step Size (delta x)");
ylabel("Maximum Error");

figure();
loglog(dx, rmserror(1, :), 'DisplayName', 'Finite Difference');
hold on;
loglog(dx, rmserror(2, :), 'DisplayName', 'Finite Element Method');
loglog(dx, rmserror(3, :), 'DisplayName', 'Spectral Method');
hold off
legend('location', 'southwest');
title("RMS Error of Approximation Methods vs Step Size");
xlabel("Step Size (delta x)");
ylabel("RMS Error");

