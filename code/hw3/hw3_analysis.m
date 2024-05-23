%% HW3 Analysis
%  Cameron Wolfe 05/08/2024

%% Loading files
x_mesh = readmatrix("problem_definition/x_mesh.mat", "FileType", "text");
y_mesh = readmatrix("problem_definition/y_mesh.mat", "FileType", "text");
boundary_mask = readmatrix("problem_definition/boundary_mask.mat", "FileType", "text");
u_boundary = complex_readmatrix("output/u2_boundary.dat");
q_boundary = complex_readmatrix("output/q2_boundary.dat");
u_sample = complex_readmatrix("output/u_sample.dat");

%% Input parsing
u_mesh = zeros(size(boundary_mask));
u_mesh(boundary_mask == 1) = u_boundary;
u_mesh(boundary_mask == 0) = u_sample;
u_mesh(boundary_mask == -1) = 1;

x_mesh = [x_mesh; x_mesh(1,:)];
y_mesh = [y_mesh; y_mesh(1,:)];
u_mesh = [u_mesh; u_mesh(1,:)];

%% Plotting setup
outer_boundary = [cos([0:pi/24:2*pi]); sin([0:pi/24:2*pi])]';
inner_boundary = 0.5 * [cos([0:pi/12:2*pi]); sin([0:pi/12:2*pi])]';

num_contours = 20;
cmap = colormap(parula(num_contours));

%% Plotting
figure(1);
colormap(cmap);
surf(x_mesh, y_mesh, real(u_mesh), 'EdgeAlpha', 0.1);
hold on;
plot3(inner_boundary(:, 1), inner_boundary(:, 2), real(u_boundary([1:24 1])) + 0.75, 'k', 'LineWidth', 2);
plot3(outer_boundary(:, 1), outer_boundary(:, 2), real(u_boundary([25:end 25])) + 0.1, 'k', 'LineWidth', 2);
hold off;
title("Helmholtz System Potential Surface");
xlabel("x");
ylabel("y");
zlabel("Potenial");
colorbar;

figure(2);
colormap(cmap);
contourf(x_mesh, y_mesh, real(u_mesh), num_contours, 'EdgeAlpha', 0.1);
hold on;
plot(inner_boundary(:, 1), inner_boundary(:, 2), 'k', 'LineWidth', 2);
plot(outer_boundary(:, 1), outer_boundary(:, 2), 'k', 'LineWidth', 2);
hold off;
axis equal;
xlabel("x");
ylabel("y");
title("Helmholtz System Potential Contour");
colorbar;

figure(3);
tiledlayout(2,1);
nexttile;
plot(0:pi/12:2*pi, real(u_boundary([1:24 1])));
hold on;
plot(0:pi/24:2*pi, real(u_boundary([25:end 25])));
hold off;
xlim([0, 2*pi])
title("Boundary Condition Values vs Theta");
ylabel("Potential");
legend(["Inner Boundary", "Outer Boundary"]);
nexttile;
plot(0:pi/12:2*pi, real(q_boundary([1:24 1])));
hold on;
plot(0:pi/24:2*pi, real(q_boundary([25:end 25])));
hold off;
ylabel("Flux");
xlabel("Theta (radians)");
legend(["Inner Boundary", "Outer Boundary"]);
xlim([0, 2*pi])

figure(4);
plot(x_mesh(1, :), real(u_mesh(1,:)));
yl = ylim;
hold on;
plot([0.5, 0.5], yl);
hold off;
title("Potential vs Radial Distance");
xlabel("r");
ylabel("Potential");
legend(["Potential", "Region Boundary"]);
