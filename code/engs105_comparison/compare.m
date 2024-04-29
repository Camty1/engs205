plotting = true;

x_fd = readmatrix("output/x080.dat");
y_fd = readmatrix("output/y080.dat");
r_fd = readmatrix("output/r080.dat");
phi_fd = readmatrix("output/phi080.dat");
u_fd = readmatrix("output/u_ss_1_080.dat");

x_fd_mesh = reshape(x_fd, 80, 80);
y_fd_mesh = reshape(y_fd, 80, 80);
u_fd_mesh = reshape(u_fd, 80, 80);
r_fd_mesh = reshape(r_fd, 80, 80);
phi_fd_mesh = reshape(phi_fd, 80, 80);
dr = r_fd_mesh(2, 1) - r_fd_mesh(1, 1);
dphi = phi_fd_mesh(1, 2) - phi_fd_mesh(1, 1);

u_fd_boundary = [u_fd_mesh(1:end-1, 1)', u_fd_mesh(end, 1:end-1), u_fd_mesh(end:-1:2, end)', u_fd_mesh(1, end:-1:2)];

dudn_fd_1 = zeros(1, 80);
dudn_fd_2 = zeros(1, 80);
dudn_fd_3 = zeros(1, 80);
dudn_fd_4 = zeros(1, 80);

dudn_fd_1 = -1 ./ r_fd_mesh(:, 1)' .* 1 / dphi .* (-3/2 * u_fd_mesh(:, 1) + 2 * u_fd_mesh(:, 2) - 1/2 * u_fd_mesh(:, 3))';
dudn_fd_2 = 1 / dr * (3/2 * u_fd_mesh(end, :) - 2 * u_fd_mesh(end-1, :) + 1/2* u_fd_mesh(end-2, :));
dudn_fd_3 = 1 ./ r_fd_mesh(end:-1:1, end)' .* 1 / dphi .* (3/2 * u_fd_mesh(end:-1:1, end) - 2 * u_fd_mesh(end:-1:1, end-1) + 1/2* u_fd_mesh(end:-1:1, end-2))';
dudn_fd_4 = -1 / dr * (-3/2 * u_fd_mesh(1, end:-1:1) + 2 * u_fd_mesh(2, end:-1:1) - 1/2 * u_fd_mesh(3, end:-1:1));

x_bem_mesh = readmatrix("x_mesh.dat");
y_bem_mesh = readmatrix("y_mesh.dat");
u_bem = readmatrix("output/u_interior.dat");

u_bem_mesh = reshape(u_bem, sqrt(length(u_bem)), []);
u_bem_boundary = readmatrix("output/u_boundary.dat");
dudn_bem_boundary = readmatrix("output/dudn_boundary.dat");


if plotting
	cmap = colormap(parula(15));
	figure(1);
	tcl = tiledlayout(1, 2);

	nexttile;
	contourf(x_fd_mesh, y_fd_mesh, u_fd_mesh, 200, 'EdgeAlpha', 0);
	colormap(cmap);
	xlim([0 1]);
	ylim([0 1]);
	clim([-1 1]);
	xlabel("x");
	ylabel("y");
	title("Finite Difference");
	colorbar;
	
	nexttile;
	contourf(x_bem_mesh, y_bem_mesh, u_bem_mesh, 200, 'EdgeAlpha', 0);
	colormap(cmap);
	xlim([0 1]);
	ylim([0 1]);
	clim([-1 1]);
	xlabel("x");
	ylabel("y");
	title("BEM");
	colorbar;

	title(tcl, "Potential Contour of Different Solution Methods");

	figure(2);
	tcl = tiledlayout(1, 2);

	nexttile;
	surf(x_fd_mesh, y_fd_mesh, u_fd_mesh, 'EdgeAlpha', 0.1);
	xlim([0 1]);
	ylim([0 1]);
	clim([-1 1]);
	xlabel("x");
	ylabel("y");
	zlabel("u");
	title("Finite Difference");
	colormap(cmap);

	nexttile;
	surf(x_bem_mesh, y_bem_mesh, u_bem_mesh, 'EdgeAlpha', 0.1);
	xlim([0 1]);
	ylim([0 1]);
	clim([-1 1]);
	xlabel("x");
	ylabel("y");
	zlabel("u");
	title("BEM");
	colormap(cmap);

	title(tcl, "Surface Plot of Different Solution Methods");

	figure(3);
	subplot(2, 1, 1);
	plot([(0:315)] * 40 / 79, u_fd_boundary);
	hold on;
	plot((1:length(u_bem_boundary))-1, u_bem_boundary);
	hold off;
	title("Boundary Potential for Finite Difference and BEM Solution");
	legend(["FD", "BEM"], "Location", "southeast");
	ylabel("Potential");

	subplot(2, 1, 2);
	plot([(0:79),(79:158), (158:237), (237:316)] * 40 / 79, [dudn_fd_1, dudn_fd_2, dudn_fd_3, dudn_fd_4]);
	hold on;
	plot((1:length(dudn_bem_boundary))-1, dudn_bem_boundary);
	hold off;
	legend(["FD", "BEM"], "Location", "southeast");
	title("Boundary Flux for Finite Difference and BEM Solution");
	xlabel("Boundary Node (FD Normalized)");
	ylabel("Flux");
end
