x_mesh = readmatrix("x.dat");
y_mesh = readmatrix("y.dat");

mask = readmatrix("mask.dat");

u_sample = complex_readmatrix("output/u_simple_sample.dat");
u = zeros(size(mask));
u(abs(mask) == 1) = 1;
u(mask == 0) = u_sample;
x_mesh(end+1, :) = x_mesh(1, :);
y_mesh(end+1, :) = y_mesh(1, :);
u(end+1, :) = u(1, :);

surf(x_mesh, y_mesh, real(u));
