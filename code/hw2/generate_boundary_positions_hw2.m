DEF_DIR = 'problem_definition';
NUM_NODE = 24;
NUM_ELEM = 24;

nodes = (0:(NUM_NODE-1))';
theta = linspace(pi, 2 * pi, 13);

x = [linspace(1, -1, 13), cos(theta(2:end-1))]';
y = [zeros(1, 13), sin(theta(2:end-1))]';
node_list = [nodes, x, y];

node_bc_type = [1, zeros(1, 11), ones(1, 12)]';
node_type_2 = 2 * sin(pi * x / 2);
node_bc_val = zeros(24, 1);
node_bc_val(node_bc_type == 2) = type_2(node_bc_type == 2);

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
writematrix([nodes, node_bc_type, node_bc_val], [DEF_DIR '/hw2.node_bcs'], 'FileType', 'text');

% New, BEMPAK like formulation
writematrix([node_list, node_bc_type, ones(24, 1), zeros(24, 2)], [DEF_DIR '/hw2.nodb'], 'FileType', 'text');
writematrix([elem_list, ones(24, 1), zeros(24, 1), elem_bc_type, elem_bc_val], [DEF_DIR '/hw2.eleb'], 'FileType', 'text');
