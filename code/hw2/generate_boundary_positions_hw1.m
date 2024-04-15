DEF_DIR = 'problem_definition';
NUM_NODE = 24;
NUM_ELEM = 24;

nodes = (0:(NUM_NODE-1))';
theta = linspace(0, 2 * pi, NUM_NODE+1);
theta = theta(1:end-1)';

x = cos(theta);
y = sin(theta);
node_list = [nodes, x, y];

node_bc_type = [1, zeros(1, 11), ones(1, 12)]';
node_type_2 = 2 * sin(pi / 2 - theta);
node_bc_val = zeros(24, 1);
node_bc_val(node_bc_type == 2) = type_2(node_bc_type == 2);

elems = (0:(NUM_ELEM-1))';
elem_list = [elems, nodes, mod(nodes+1, NUM_NODE)];

elem_thetas = theta + theta(1);

elem_bc_type = [2 * ones(1, 12), zeros(1, 12)]';
elem_bc_type_2 = 2 * sin(pi / 2 - elem_thetas);
elem_bc_val = zeros(24, 1);
elem_bc_val(elem_bc_type == 2) = elem_bc_type_2(elem_bc_type == 2);

% For old formulation
writematrix(node_list ,[DEF_DIR '/hw1.nod'], 'FileType', 'text');
writematrix(elem_list, [DEF_DIR '/hw1.ele'], 'FileType', 'text');
writematrix([nodes, node_bc_type, node_bc_val], [DEF_DIR '/hw1.node_bcs'], 'FileType', 'text');

% New, BEMPAK like formulation
writematrix([node_list, node_bc_type, ones(24, 1), zeros(24, 2)], [DEF_DIR '/hw1.nodb'], 'FileType', 'text');
writematrix([elem_list, ones(24, 1), zeros(24, 1), elem_bc_type, elem_bc_val], [DEF_DIR '/hw1.eleb'], 'FileType', 'text');
