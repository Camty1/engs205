r_sample = 0.025:0.025:0.75;
theta = 0:pi/12:23*pi/12;

[r_sample_mesh, theta_mesh] = meshgrid(r_sample, theta);

x_sample_mesh = r_sample_mesh .* cos(theta_mesh);
y_sample_mesh = r_sample_mesh .* sin(theta_mesh);

r_sample_out = [0; reshape(r_sample_mesh, [], 1)];
theta_out = [0; reshape(theta_mesh, [], 1)];

sample_positions = [r_sample_out .* cos(theta_out), r_sample_out .* sin(theta_out)];

writematrix([(0:length(r_sample_out)-1)' sample_positions], "problem_definition/domain_samples.dat");
