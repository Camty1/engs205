num_r = 41;
num_theta = 48;
num_node = 24;


bc_theta = 0:2*pi/num_node:2*pi*(1 - 1 / num_node);
r = linspace(0, 1, num_r);
theta = 0:2*pi/num_theta:2*pi*(1 - 1/num_theta);

[r_mesh, theta_mesh] = meshgrid(r, theta);

mask = zeros(size(r_mesh));

mask(r_mesh == 1 & ismember(round(theta_mesh * num_theta / pi), round(bc_theta * num_theta / pi))) = 1;
mask(r_mesh == 1 & ~ismember(round(theta_mesh * num_theta / pi), round(bc_theta * num_theta / pi))) = -1; % outside boundary

x_mesh = r_mesh .* cos(theta_mesh);
y_mesh = r_mesh .* sin(theta_mesh);

x_boundary = x_mesh(mask == 1)';
y_boundary = y_mesh(mask == 1)';
x_sample = x_mesh(mask == 0)';
y_sample = y_mesh(mask == 0)';

num_sample = length(x_sample);

writematrix([0:(num_node-1); x_boundary; y_boundary]', "problem_definition/simple.nod", 'FileType', 'text');
writematrix([(0:num_node-1); (0:num_node-1); mod(1:num_node, num_node); ones(1, num_node); zeros(1, num_node)]', "problem_definition/simple.ele", 'FileType', 'text');
writematrix([0:(num_node-1); ones(1, num_node); ones(1, num_node)]', "problem_definition/simple.bcs", 'FileType', 'text');
writematrix([0:(num_sample-1); x_sample; y_sample; ones(1,num_sample)]', "problem_definition/simple.spl", 'FileType', 'text');

writematrix(x_mesh, "x.dat");
writematrix(y_mesh, "y.dat");
writematrix(mask, "mask.dat");
