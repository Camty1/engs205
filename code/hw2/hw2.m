NUM_NODE = 24;
node_pos = readmatrix("problem_definition/hw2.nod", 'FileType', 'text');
x_boundary = node_pos(:, 2);
y_boundary = node_pos(:, 3);

u_boundary = readmatrix("output/u_boundary_hw2c.dat");
dudn_boundary = readmatrix("output/dudn_boundary_hw2c.dat");

interior_pos = readmatrix("problem_definition/sample_points_hw2.nod", 'FileType', 'text');
x_interior = interior_pos(:, 2);
y_interior = interior_pos(:, 3);

right_corner_x = (0:0.01:1)';
right_corner_y = (-0.5:0.005:0)';
left_corner_x = -right_corner_x;
left_corner_y = (-0.5:0.005:0)';

u_interior = readmatrix("output/u_interior_hw2c.dat");
dudx_interior = readmatrix("output/dudx_interior_hw2c.dat");
dudy_interior = readmatrix("output/dudy_interior_hw2c.dat");

x_mesh = readmatrix("x_mesh_hw2.mat", 'FileType', 'text');
y_mesh = readmatrix("y_mesh_hw2.mat", 'FileType', 'text');
mask = readmatrix("mask_hw2.mat", 'FileType', 'text');

u_mesh = nan(size(x_mesh));
dudx_mesh = nan(size(x_mesh));
dudy_mesh = nan(size(x_mesh));

boundary = [[1; cos(pi:pi/24:2*pi)'], [0; sin(pi:pi/24:2*pi)']];

u_mesh(mask == 1) = u_interior;
dudx_mesh(mask == 1) = dudx_interior;
dudy_mesh(mask == 1) = dudy_interior;

u_right_corner = interp2(x_mesh, y_mesh, u_mesh, right_corner_x, right_corner_y);
u_left_corner = interp2(x_mesh, y_mesh, u_mesh, left_corner_x, left_corner_y);

close all;
figure(1);
plot(u_boundary);
hold on;
plot(dudn_boundary);
hold off;
title("Boundary Potential and Flux by Node");
xlabel("Node (n)");
ylabel("Magnitude");
legend("Boundary Potential", "Boundary Flux", "location", "northeast");
xlim([1, 24]);

figure(2);
contourf(x_mesh, y_mesh, u_mesh);
hold on;
plot(boundary(:, 1), boundary(:, 2), 'k', 'LineWidth', 2.5)
hold off;
axis equal
colorbar
title("Contour of Potential of Interior");
xlabel("x");
ylabel("y");


figure(3);
quiver(x_mesh, y_mesh, dudx_mesh, dudy_mesh, 3, 'LineWidth', 2);
hold on;
plot(boundary(:, 1), boundary(:, 2), 'k', 'LineWidth', 2.5)
hold off;
ylim([-1.125, 0.125])
axis equal;
title("Gradient of Potential");
xlabel("x")
ylabel("y")

figure(4);
plot(sqrt((right_corner_x - 1).^2 + right_corner_y.^2), u_right_corner);
hold on;
plot(sqrt((left_corner_x + 1).^2 + left_corner_y.^2), u_left_corner);
hold off;
title("Potential vs Distance from Corner");
xlabel("Distance");
ylabel("Potential");
legend("Right Corner", "Left Corner");
