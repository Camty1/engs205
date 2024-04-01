NUM_NODE = 24;
node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', 'text');
bcs = readmatrix("problem_definition/hw1.bcs", 'FileType', 'text');
new_bcs = readmatrix("output/new_bcs.dat");
u = zeros(NUM_NODE, 1);
dudn = zeros(NUM_NODE, 1);

for i=1:NUM_NODE
	if bcs(i, 2) == 1
		u(i) = bcs(i, 3);
		dudn(i) = new_bcs(i);
	else
		dudn(i) = bcs(i, 3);
		u(i) = new_bcs(i);
	end
end

close all;
figure(1);
plot(1:NUM_NODE, u);
hold on;
plot(1:NUM_NODE, dudn);
hold off;
