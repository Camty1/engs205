NUM_ROWS = 31;
NUM_COLS = 31;

x = linspace(-1, 1, NUM_COLS);
y = linspace(-1, 0, NUM_ROWS);
x_shift = (x(2:end) + x(1:end-1)) / 2;
y_shift = (y(2:end) + y(1:end-1)) / 2;

node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', "text");
elem_list = readmatrix("problem_definition/hw1.ele", 'FileType', "text");
x_boundary = node_pos(:, 2);
y_boundary = node_pos(:, 3);

midpoint_pos = [(x_boundary(elem_list(:, 2) + 1) + x_boundary(elem_list(:, 3) + 1)) / 2, (y_boundary(elem_list(:, 2) + 1) + y_boundary(elem_list(:, 3) + 1)) / 2];

[x_mesh, y_mesh] = meshgrid(x_shift,y_shift);

r_mesh = x_mesh.^2 + y_mesh.^2;

mask = ones(size(x_mesh));

for i=1:size(x_mesh, 1)
	for j=1:size(x_mesh, 2)
		dot_best = 0;
		l_best = -1;
		r = [x_mesh(i, j), y_mesh(i, j)];
		for l=1:24
			dot_val = sum(r .* midpoint_pos(l, :));
			if dot_val > dot_best
				dot_best = dot_val;
				l_best = l;
			end
		end
		r_minus = r - midpoint_pos(l_best, :);
		if sum(r_minus .* midpoint_pos(l_best, :)) >= 0
			mask(i, j) = 0;
		end
	end
end


mask_2 = ones(size(x_mesh));

mask_2(r_mesh >= 0.95) = 0;
mask_2(y_mesh >= -0.025) = 0;

mask = mask & mask_2;

x_out = [reshape(x_mesh(mask ~= 0), [], 1)];
y_out = [reshape(y_mesh(mask ~= 0), [], 1)];

writematrix([(0:length(x_out)-1)', x_out, y_out], "problem_definition/sample_points_hw2.nod", 'FileType', 'text');

writematrix(x_mesh, "x_mesh_hw2.mat", 'FileType', 'text');
writematrix(y_mesh, "y_mesh_hw2.mat", 'FileType', 'text');
writematrix(mask, "mask_hw2.mat", 'FileType', 'text');
