elem_list = readmatrix("problem_definition/hw1.ele", 'FileType', 'text');
node_pos = readmatrix("problem_definition/hw1.nod", 'FileType', 'text');

for i=1:24

	plot(node_pos(elem_list(i, 2:3) + 1, 2), node_pos(elem_list(i, 2:3) + 1, 3) )
	if i == 1
		hold on;
	end
end
