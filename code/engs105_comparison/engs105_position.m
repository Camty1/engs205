r_min = 0.1;
r_max = 1;
theta_min = 0;
theta_max = pi / 2;

nodes_per_side = 41;
sample_points_per_axis = 101;

r = linspace(r_min, r_max, nodes_per_side);
theta = linspace(theta_min, theta_max, nodes_per_side);

x = [r(1:end-1), r_max * cos(theta(1:end-1)), zeros(1, nodes_per_side-1), r_min * cos(theta(end:-1:2))]';
y = [zeros(1, nodes_per_side-1), r_max * sin(theta(1:end-1)), r(end:-1:2), r_min * sin(theta(end:-1:2))]';

nodes = [(1:length(x))'-1, x, y];
elems = [(1:length(x))'-1, (1:length(x))'-1, mod(1:length(x), length(x))'];

bc_types = [2 * ones(1, nodes_per_side - 1), ones(1, nodes_per_side-1), ones(1, nodes_per_side), 2*ones(1, nodes_per_side-2)]';

bc_vals = [ones(1,nodes_per_side-1), - r_max * cos(4*theta(1:end-1)), -ones(1,nodes_per_side), ones(1, nodes_per_side-2)]';

bcs = [nodes(:, 1), bc_types, bc_vals];

writematrix(nodes, "problem_definition/hw1.nod", 'FileType', 'text');
writematrix(elems, "problem_definition/hw1.ele", 'FileType', 'text');
writematrix(bcs, "problem_definition/hw1.bcs", 'FileType', 'text');

r_sample = linspace(r_min, r_max, sample_points_per_axis);
theta_sample = linspace(theta_min, theta_max, sample_points_per_axis);

[r_mesh, theta_mesh] = meshgrid(r_sample, theta_sample);

x_mesh = r_mesh .* cos(theta_mesh);
y_mesh = r_mesh .* sin(theta_mesh);

exclusion = 0.01;

x_mesh = x_mesh * (1 - 2 * exclusion) + exclusion;
y_mesh = y_mesh * (1 - 2 * exclusion) + exclusion;

x_sample = reshape(x_mesh, [], 1);
y_sample = reshape(y_mesh, [], 1);

sample_points = [(1:length(x_sample))'-1, x_sample, y_sample];

writematrix(sample_points, "problem_definition/sample_points.nod", 'FileType', 'text');

writematrix(x_mesh, "x_mesh.dat");
writematrix(y_mesh, "y_mesh.dat");

plot(nodes(bc_types == 1, 2), nodes(bc_types == 1, 3), 'x');
hold on;
plot(nodes(bc_types == 2, 2), nodes(bc_types == 2, 3), 'x');
plot(x_mesh, y_mesh, '.')
hold off;


