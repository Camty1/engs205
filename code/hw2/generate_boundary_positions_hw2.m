DEF_DIR = 'problem_definition';
NUM_NODE = 24;
NUM_ELEM = 24;

nodes = (0:(NUM_NODE-1))';
theta = linspace(pi, 2 * pi, 13);

elems = (0:(NUM_ELEM-1))';
elem_list = [elems, nodes, mod(nodes+1, NUM_NODE)];
elem_bc_type = [2 * ones(1, 12), ones(1, 12)]';


x = [linspace(1, -1, 13), cos(theta(2:end-1))]';
y = [zeros(1, 13), sin(theta(2:end-1))]';
node_list = [nodes, x, y];

node_bc_type_old = [1; 2 * ones(11,1); ones(12, 1)];
node_bc_type_new = zeros(NUM_NODE, 1);
for i=1:NUM_ELEM
	node_bc_type_new(elem_list(i, 2) + 1) = node_bc_type_new(elem_list(i, 2) + 1) + elem_bc_type(i);
	node_bc_type_new(elem_list(i, 3) + 1) = 10 * elem_bc_type(i) + node_bc_type_new(elem_list(i, 3) + 1);
end
node_type_2 = 2 * sin(pi * x / 2);
node_bc_val_old = zeros(24, 1);
node_bc_val_old(node_bc_type_old == 2) = node_type_2(node_bc_type_old == 2);
node_bc_val_new = zeros(24, 2);
node_bc_val_new(node_bc_type_new == 22, :) = node_type_2(node_bc_type_new == 22) * ones(1,2);
node_bc_val_new(node_bc_type_new == 21, :) = [-2, 0];
node_bc_val_new(node_bc_type_new == 12, :) = [0, 2];
node_bc_val_new(node_bc_type_new == 11, :) = ones(sum(node_bc_type_new == 11), 1) * [0, nan(1,1)];


elems = (0:(NUM_ELEM-1))';
elem_list = [elems, nodes, mod(nodes+1, NUM_NODE)];

elem_midpoints = [(x(elem_list(:, 2)+1) + x(elem_list(:, 3)+1))/2, (y(elem_list(:, 2)+1) + y(elem_list(:, 3)+1))/2];

elem_bc_type = [2 * ones(1, 12), zeros(1, 12)]';
elem_bc_type_2 = 2 * sin(pi * elem_midpoints(:, 1) / 2);
elem_bc_val = zeros(24, 1);
elem_bc_val(elem_bc_type == 2) = elem_bc_type_2(elem_bc_type == 2);

% For old formulation
writematrix(node_list ,[DEF_DIR '/hw2.nod'], 'FileType', 'text');
writematrix(elem_list, [DEF_DIR '/hw2.ele'], 'FileType', 'text');
writematrix([nodes, node_bc_type_old, node_bc_val_old], [DEF_DIR '/hw2.node_bcs'], 'FileType', 'text');

% New, BEMPAK like formulation
writematrix([node_list, node_bc_type_new, node_bc_val_new], [DEF_DIR '/hw2.nodb'], 'FileType', 'text');
writematrix([elem_list, ones(24, 1), zeros(24, 1), elem_bc_type, elem_bc_val], [DEF_DIR '/hw2.eleb'], 'FileType', 'text');
