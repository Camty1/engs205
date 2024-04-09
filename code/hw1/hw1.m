NUM_NODE = 24;
node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', 'text');

u_boundary = readmatrix("u.dat");
dudn_boundary = readmatrix("dudn.dat");

close all;
figure(1);
plot(u_boundary);
hold on;
plot(dudn_boundary);
hold off;
