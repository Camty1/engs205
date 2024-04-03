NUM_NODE = 24;
node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', 'text');

u_boundary = readmatrix("output/u_boundary.dat");
dudn_boundary = readmatrix("output/dudn_boundary.dat");

close all;
figure(1);
plot(u_boundary);
hold on;
plot(dudn_boundary);
hold off;
