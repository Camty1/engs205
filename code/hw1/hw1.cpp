#include "hw1.hpp"
#include "utils_205.hpp"
#define NUM_NODE 24
#define NUM_ELEM 24
#define QUAD_DEG 5
#define NUM_SRC 1
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433

int main(int argc, char **argv) {
  // Matrix and vector definitions
  double *node_pos = new double[2 * NUM_NODE];

  int *elem_list = new int[2 * NUM_ELEM];

  int *bc_types = new int[NUM_NODE];
  double *bc_values = new double[NUM_NODE];

  double *src_pos = new double[2 * NUM_SRC];
  double *src_values = new double[NUM_SRC];

  double *A = new double[NUM_NODE * NUM_NODE];
  fill('c', A, NUM_NODE, NUM_NODE, 0, NUM_NODE - 1, 0, NUM_NODE - 1, 0.0);
  double *B = new double[NUM_NODE * NUM_NODE];
  fill('c', B, NUM_NODE, NUM_NODE, 0, NUM_NODE - 1, 0, NUM_NODE - 1, 0.0);
  double *F = new double[NUM_NODE];
  fill(F, NUM_NODE, 0, NUM_NODE - 1, 0.0);

  double *quad_points = new double[QUAD_DEG];
  double *quad_weights = new double[QUAD_DEG];

  double *x_elem = new double[2 * NUM_ELEM];
  double *y_elem = new double[2 * NUM_ELEM];
  double *delta_s = new double[NUM_ELEM];

  double alpha = 11 * PI / 24.0;

  // Filename definitions
  std::string pos_filename = "problem_definition/hw1.nod";
  std::string elem_filename = "problem_definition/hw1.ele";
  std::string bc_filename = "problem_definition/hw1.bcs";
  std::string src_filename = "problem_definition/hw1.src";
  std::string quad_filename =
      "quadrature/deg_" + std::to_string(QUAD_DEG) + ".gqd";

  // Read from files
  // Nodes
  std::ifstream pos_file(pos_filename);

  if (!pos_file.is_open()) {
    std::cerr << "Could not open position file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_NODE; i++) {
    std::string line;
    std::string token;
    getline(pos_file, line);
    std::istringstream lineStream(line);

    // Node number, do nothing
    getline(lineStream, token, ',');

    // Node x
    getline(lineStream, token, ',');
    double x = std::stod(token);

    int idx = general_matrix_index('c', NUM_NODE, 2, i, 0);
    node_pos[idx] = x;

    // Node y
    getline(lineStream, token, ',');
    double y = std::stod(token);

    idx = general_matrix_index('c', NUM_NODE, 2, i, 1);
    node_pos[idx] = y;
  }

  // Elements
  std::ifstream elem_file(elem_filename);

  if (!elem_file.is_open()) {
    std::cerr << "Could not open element file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_ELEM; i++) {
    std::string line;
    std::string token;
    getline(elem_file, line);
    std::istringstream lineStream(line);

    // Element number, do nothing
    getline(lineStream, token, ',');

    // First node
    getline(lineStream, token, ',');
    int idx = general_matrix_index('c', NUM_ELEM, 2, i, 0);
    elem_list[idx] = std::stoi(token);

    // Second node
    getline(lineStream, token, ',');
    idx = general_matrix_index('c', NUM_ELEM, 2, i, 1);
    elem_list[idx] = std::stoi(token);
  }

  // BCs
  std::ifstream bc_file(bc_filename);

  if (!bc_file.is_open()) {
    std::cerr << "Could not open bc file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_NODE; i++) {
    std::string line;
    std::string token;
    getline(bc_file, line);
    std::istringstream lineStream(line);

    // Node number, do nothing
    getline(lineStream, token, ',');

    // BC type
    getline(lineStream, token, ',');
    bc_types[i] = std::stoi(token);

    // BC value
    getline(lineStream, token, ',');
    bc_values[i] = std::stod(token);
  }

  // Sources
  std::ifstream src_file(src_filename);

  if (!src_file.is_open()) {
    std::cerr << "Could not open src file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_SRC; i++) {
    std::string line;
    std::string token;
    getline(src_file, line);
    std::istringstream lineStream(line);

    // Source number, do nothing
    getline(lineStream, token, ',');

    // Source x pos
    getline(lineStream, token, ',');
    int idx = general_matrix_index('c', NUM_SRC, 2, i, 0);
    src_pos[idx] = std::stod(token);

    // Source y pos
    getline(lineStream, token, ',');
    idx = general_matrix_index('c', NUM_SRC, 2, i, 1);
    src_pos[idx] = std::stod(token);

    // Source value
    getline(lineStream, token, ',');
    src_values[i] = std::stod(token);
  }

  // Quadrature
  std::ifstream quad_file(quad_filename);

  if (!quad_file.is_open()) {
    std::cerr << "Could not open quadrature file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < QUAD_DEG; i++) {
    std::string line;
    std::string token;
    getline(quad_file, line);
    std::istringstream lineStream(line);

    // Point number, do nothing
    getline(lineStream, token, ',');

    // Quadrature point
    getline(lineStream, token, ',');
    quad_points[i] = std::stod(token);

    // Quadrature weight
    getline(lineStream, token, ',');
    quad_weights[i] = std::stod(token);
  }

  // Print inputs
  if (PRINT_INPUTS) {
    std::cout << "Node positions: " << std::endl;
    print_matrix('c', node_pos, NUM_NODE, 2);

    std::cout << "Element list: " << std::endl;
    print_matrix('c', elem_list, NUM_ELEM, 2);

    std::cout << "BC types: " << std::endl;
    print_vector(bc_types, NUM_NODE);

    std::cout << "BC values: " << std::endl;
    print_vector(bc_values, NUM_NODE);

    std::cout << "Source positions: " << std::endl;
    print_matrix('c', src_pos, NUM_SRC, 2);

    std::cout << "Source values: " << std::endl;
    print_vector(src_values, NUM_SRC);

    std::cout << "Quadrature points: " << std::endl;
    print_vector(quad_points, QUAD_DEG);

    std::cout << "Quadrature weights: " << std::endl;
    print_vector(quad_weights, QUAD_DEG);
  }

  // Get element node position matrices
  for (int i = 0; i < NUM_ELEM; i++) {
    for (int j = 0; j < 2; j++) {
      int elem_idx = general_matrix_index('c', NUM_ELEM, 2, i, j);
      int x_idx =
          general_matrix_index('c', NUM_NODE, 2, elem_list[elem_idx], 0);
      int y_idx =
          general_matrix_index('c', NUM_NODE, 2, elem_list[elem_idx], 1);
      x_elem[elem_idx] = node_pos[x_idx];
      y_elem[elem_idx] = node_pos[y_idx];
    }
    delta_s[i] = calc_delta_s(x_elem, y_elem, NUM_ELEM, i);
  }

  // Populating A and B
  for (int i = 0; i < NUM_NODE; i++) {
    for (int j = 0; j < NUM_ELEM; j++) {
      int elem_idx_1 = general_matrix_index('c', NUM_ELEM, 2, j, 0);
      int elem_idx_2 = general_matrix_index('c', NUM_ELEM, 2, j, 1);

      if (i == elem_list[elem_idx_1] || i == elem_list[elem_idx_2]) {
        int mat_index_1_1 =
            general_matrix_index('c', NUM_NODE, NUM_NODE, elem_list[elem_idx_1],
                                 elem_list[elem_idx_1]);
        int mat_index_2_1 =
            general_matrix_index('c', NUM_NODE, NUM_NODE, elem_list[elem_idx_2],
                                 elem_list[elem_idx_1]);
        int mat_index_1_2 =
            general_matrix_index('c', NUM_NODE, NUM_NODE, elem_list[elem_idx_1],
                                 elem_list[elem_idx_2]);
        int mat_index_2_2 =
            general_matrix_index('c', NUM_NODE, NUM_NODE, elem_list[elem_idx_2],
                                 elem_list[elem_idx_2]);

        A[mat_index_1_1] += alpha / 2.0;
        A[mat_index_2_2] += alpha / 2.0;

        B[mat_index_1_1] += delta_s[j] * (1.5 - log(delta_s[j])) / 2.0;
        B[mat_index_1_2] += delta_s[j] * (0.5 - log(delta_s[j])) / 2.0;
        B[mat_index_2_1] += delta_s[j] * (0.5 - log(delta_s[j])) / 2.0;
        B[mat_index_2_2] += delta_s[j] * (1.5 - log(delta_s[j])) / 2.0;

      } else {
        for (int k = 0; k < QUAD_DEG; k++) {
          double eta = quad_points[k];
          double w = quad_weights[k];
          double *phi = new double[2];

          calc_basis(eta, phi);

          double x_elem_1 = x_elem[elem_idx_1];
          double x_elem_2 = x_elem[elem_idx_2];
          double y_elem_1 = y_elem[elem_idx_1];
          double y_elem_2 = y_elem[elem_idx_2];

          double xs = x_elem_1 * phi[0] + x_elem_2 * phi[1];
          double ys = y_elem_1 * phi[0] + y_elem_2 * phi[1];

          int node_x_idx = general_matrix_index('c', NUM_NODE, 2, i, 0);
          int node_y_idx = general_matrix_index('c', NUM_NODE, 2, i, 1);

          double x_node = node_pos[node_x_idx];
          double y_node = node_pos[node_y_idx];

          double r = sqrt(pow(xs - x_node, 2) + pow(ys - y_node, 2));

          double drdn = ((y_elem_2 - y_elem_1) * (xs - x_node) -
                         (x_elem_2 - x_elem_1) * (ys - y_node)) /
                        (delta_s[j] * r);

          double G = -1 / (2 * PI) * log(r);
          double dGdn = -1 / (2 * PI * r) * drdn;

          int mat_idx_1 = general_matrix_index('c', NUM_NODE, NUM_NODE, i,
                                               elem_list[elem_idx_1]);
          int mat_idx_2 = general_matrix_index('c', NUM_NODE, NUM_NODE, i,
                                               elem_list[elem_idx_2]);
          A[mat_idx_1] += phi[0] * dGdn * delta_s[j] * w / 2.0;
          A[mat_idx_2] += phi[1] * dGdn * delta_s[j] * w / 2.0;

          B[mat_idx_1] += phi[0] * G * delta_s[j] * w / 2.0;
          B[mat_idx_2] += phi[1] * G * delta_s[j] * w / 2.0;
        }
      }
    }
  }

  // Populate F
  for (int i = 0; i < NUM_NODE; i++) {
    int node_x_idx = general_matrix_index('c', NUM_NODE, 2, i, 0);
    int node_y_idx = general_matrix_index('c', NUM_NODE, 2, i, 1);

    double node_x = node_pos[node_x_idx];
    double node_y = node_pos[node_y_idx];

    for (int s = 0; s < NUM_SRC; s++) {
      int src_x_idx = general_matrix_index('c', NUM_SRC, 2, s, 0);
      int src_y_idx = general_matrix_index('c', NUM_SRC, 2, s, 1);

      double src_x = src_pos[src_x_idx];
      double src_y = src_pos[src_y_idx];

      double r = sqrt(pow(src_x - node_x, 2) + pow(src_y - node_y, 2));

      F[i] -= src_values[i] / (2 * PI) * log(r);
    }
  }
  write_matrix('c', A, NUM_NODE, NUM_NODE, "A.dat");
  write_matrix('c', B, NUM_NODE, NUM_NODE, "B.dat");
  write_vector(F, NUM_NODE, "F.dat");

  return EXIT_SUCCESS;
}

double calc_delta_s(double *x_elem, double *y_elem, int num_rows, int elem) {
  int idx_1 = general_matrix_index('c', num_rows, 2, elem, 0);
  int idx_2 = general_matrix_index('c', num_rows, 2, elem, 1);

  return sqrt(pow(x_elem[idx_2] - x_elem[idx_1], 2) +
              pow(y_elem[idx_2] - y_elem[idx_1], 2));
}

int calc_basis(double eta, double *phi) {
  phi[0] = (1 - eta) / 2.0;
  phi[1] = (1 + eta) / 2.0;

  return EXIT_SUCCESS;
}
