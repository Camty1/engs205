%% HW3 Setup
%  Cameron Wolfe 05/08/24

%% Parameters
NUM_INNER_NODES = 24;
NUM_OUTER_NODES = 48;
R_SPACING = 0.025;
THETA_SPACING = pi / 48;

r1 = 0.5;
r2 = 1;
theta_1 = linspace(0, 2 * pi, NUM_INNER_NODES + 1);
theta_2 = linspace(0, 2 * pi, NUM_OUTER_NODES + 1);

%% Sample and boundary mesh definitions
r_sample = 0:R_SPACING:r2;
theta_sample = 0:THETA_SPACING:2*pi;
theta_sample = theta_sample(1:end-1);

[r_mesh, theta_mesh] = meshgrid(r_sample, theta_sample);

x_mesh = r_mesh .* cos(theta_mesh);
y_mesh = r_mesh .* sin(theta_mesh);

boundary_mask = zeros(size(r_mesh));
boundary_mask(ismember(round(theta_sample / THETA_SPACING), round(theta_1 / THETA_SPACING)), r_sample == r1) = 1;
boundary_mask(ismember(round(theta_sample / THETA_SPACING), round(theta_2 / THETA_SPACING)), r_sample == r2) = 1;
boundary_mask(boundary_mask(:, end) == 0, end) = -1;

region_mask = ones(size(r_mesh));
region_mask(r_mesh >= r1) = 2;

x_boundary = x_mesh(boundary_mask == 1);
y_boundary = y_mesh(boundary_mask == 1);

x_sample = x_mesh(boundary_mask == 0);
y_sample = y_mesh(boundary_mask == 0);
region_sample = region_mask(boundary_mask == 0);

%% Elem list creation
inner_elems = [(0:NUM_INNER_NODES-1)', (0:NUM_INNER_NODES-1)', mod((1:NUM_INNER_NODES)', NUM_INNER_NODES)];

outer_elems = [(0:NUM_OUTER_NODES-1)', (0:NUM_OUTER_NODES-1)', mod((1:NUM_OUTER_NODES)', NUM_OUTER_NODES)] + NUM_INNER_NODES;

left = [ones(NUM_INNER_NODES, 1); 2*ones(NUM_OUTER_NODES, 1)];
right = [2*ones(NUM_INNER_NODES, 1); zeros(NUM_OUTER_NODES, 1)];

%% BC list creation
bc_type = [zeros(NUM_INNER_NODES, 1); ones(NUM_OUTER_NODES, 1)];
bc_val = [nan(NUM_INNER_NODES, 1); ones(NUM_OUTER_NODES, 1)];


%% Output to files
writematrix(x_mesh, "problem_definition/x_mesh.mat", 'FileType', 'text');
writematrix(y_mesh, "problem_definition/y_mesh.mat", 'FileType', 'text');
writematrix(boundary_mask, "problem_definition/boundary_mask.mat", 'FileType', 'text');
writematrix([((1:NUM_INNER_NODES+NUM_OUTER_NODES)-1)', x_boundary, y_boundary], "problem_definition/hw3.nod", 'FileType', 'text');
writematrix([((1:length(x_sample))-1)', x_sample, y_sample, region_sample], "problem_definition/hw3.spl", 'FileType', 'text');
writematrix([[inner_elems; outer_elems] left right], "problem_definition/hw3.ele", 'FileType', 'text');
writematrix([(0:NUM_INNER_NODES+NUM_OUTER_NODES-1)' bc_type bc_val], "problem_definition/hw3.bcs", 'FileType', 'text');

