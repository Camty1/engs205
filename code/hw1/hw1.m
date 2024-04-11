NUM_NODE = 24;
node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', 'text');
x_boundary = node_pos(:, 2);
y_boundary = node_pos(:, 3);

u_boundary = readmatrix("output/u_boundary.dat");
dudn_boundary = readmatrix("output/dudn_boundary.dat");

interior_pos = readmatrix("problem_definition/sample_points.nod", 'FileType', 'text');
x_interior = interior_pos(:, 2);
y_interior = interior_pos(:, 3);

u_interior = readmatrix("output/u_interior.dat");
dudx_interior = readmatrix("output/dudx_interior.dat");
dudy_interior = readmatrix("output/dudy_interior.dat");

x_mesh = readmatrix("x_mesh.mat", 'FileType', 'text');
y_mesh = readmatrix("y_mesh.mat", 'FileType', 'text');
mask = readmatrix("mask.mat", 'FileType', 'text');

u_mesh = nan(size(x_mesh));
dudx_mesh = nan(size(x_mesh));
dudy_mesh = nan(size(x_mesh));

u_mesh(mask == 1) = u_interior;
dudx_mesh(mask == 1) = dudx_interior;
dudy_mesh(mask == 1) = dudy_interior;

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
contourf(x_mesh, y_mesh, u_mesh, -1:0.2:1);
axis equal
colorbar
title("Contour of Potential of Interior");
xlabel("x");
ylabel("y");

figure(3);
quiver(x_mesh, y_mesh, dudx_mesh, dudy_mesh, 10);
axis equal
