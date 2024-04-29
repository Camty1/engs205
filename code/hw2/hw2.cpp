#include "hw2.hpp"
#define NUM_NODE 24
#define NUM_ELEM 24
#define NUM_QUAD 4
#define NUM_SRC 1
#define NUM_SAMP 300
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433

int main(int argc, char** argv) {
  // Define matrices and vectors
  Array<double> node_pos(2, NUM_NODE);
  Array<int> node_bc_type(NUM_NODE);
  Array<double> node_bc_vals(2, NUM_NODE);
  
  Array<int> elem_list(2, NUM_ELEM);
  Array<int> elem_boundaries(2, NUM_ELEM);

  Array<double> source_pos(2, NUM_SRC);
  Array<double> source_values(NUM_SRC);

  Array<double> quad_points(NUM_QUAD);
  Array<double> quad_weights(NUM_QUAD);

  Array<double> x_elem(2, NUM_ELEM);
  Array<double> y_elem(2, NUM_ELEM);

  Array<double> elem_lengths(NUM_ELEM);
  Array<double> elem_tangent(2, NUM_ELEM);
  Array<double> elem_normal(2, NUM_ELEM);

  Array<double> node_angles(NUM_NODE);

  Array<double> phi(2, NUM_QUAD);
  Array<double> xk(NUM_QUAD);
  Array<double> yk(NUM_QUAD);

  Array<double> G(NUM_QUAD);
  Array<double> dGdn(NUM_QUAD);
  
  Array<double> A(NUM_NODE, NUM_NODE);
  Array<double> B(NUM_NODE, 2 * NUM_NODE);
  Array<double> LHS(NUM_NODE, NUM_NODE);
  A.fill(0.0);
  B.fill(0.0);
  LHS.fill(0.0);

  Array<double> RHS(NUM_NODE);
  RHS.fill(0.0);

  // Filenames
  std::string node_filename = "problem_definition/hw2.nodb";
  std::string elem_filename = "problem_definition/hw2.ele";
  std::string source_filename = "problem_definition/hw2.src";
  std::string quad_filename =
      "quadrature/deg_" + std::to_string(NUM_QUAD) + ".gqd";

  // Read files
  // Nodes
  std::ifstream node_file(node_filename);

  if (!node_file.is_open()) {
    throw std::invalid_argument("Node file (" + node_filename + ") could not be opened, check to make sure its spelled correctly");
  }

  for (int i = 0; i < NUM_NODE; i++) {
    std::string line;
    std::string token;
    getline(node_file, line);
    std::istringstream lineStream(line);

    // Node number, do nothing
    getline(lineStream, token, ',');

    // Node x and y
    for (int j = 0; j < 2; j++) {
      getline(lineStream, token, ',');
      node_pos(j, i) = std::stod(token);
    }

    // BC type
    getline(lineStream, token, ',');
    node_bc_type(i) = std::stoi(token);

    // Node BC coeffs
    for (int j = 0; j < 2; j++) {
      getline(lineStream, token, ',');
      node_bc_vals(j, i) = std::stod(token);
    }
  }

  node_file.close();

  // Elements
  std::ifstream elem_file(elem_filename);

  if (!elem_file.is_open()) {
    throw std::invalid_argument("elem file (" + elem_filename + ") could not be opened, check to make sure its spelled correctly");
  }

  for (int i = 0; i < NUM_ELEM; i++) {
    std::string line;
    std::string token;
    getline(elem_file, line);
    std::istringstream lineStream(line);

    // Element number, do nothing
    getline(lineStream, token, ',');

    // Element incidence list
    for (int j = 0; j < 2; j++) {
      getline(lineStream, token, ',');
      elem_list(j, i) = std::stoi(token);
    }

  }

  elem_file.close();

  // Sources
  std::ifstream source_file(source_filename);

  if (!source_file.is_open()) {
    std::cerr << "Could not open source file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_SRC; i++) {
    std::string line;
    std::string token;
    getline(source_file, line);
    std::istringstream lineStream(line);

    // Source number, do nothing
    getline(lineStream, token, ',');

    // Source x pos
    getline(lineStream, token, ',');
    source_pos(0, i) = std::stod(token);

    // Source y pos
    getline(lineStream, token, ',');
    source_pos(1, i) = std::stod(token);

    // Source value
    getline(lineStream, token, ',');
    source_values(i) = std::stod(token);
  }

  source_file.close();

  // Quadrature
  std::ifstream quad_file(quad_filename);

  if (!quad_file.is_open()) {
    std::cerr << "Could not open quadrature file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_QUAD; i++) {
    std::string line;
    std::string token;
    getline(quad_file, line);
    std::istringstream lineStream(line);

    // Point number, do nothing
    getline(lineStream, token, ',');

    // Quadrature point
    getline(lineStream, token, ',');
    quad_points(i) = std::stod(token);

    // Quadrature weight
    getline(lineStream, token, ',');
    quad_weights(i) = std::stod(token);
  }

  quad_file.close();
  if (PRINT_INPUTS) {
    std::cout << "node_pos" << std::endl;
    node_pos.print();
    std::cout << "elem_list" << std::endl;;
    elem_list.print();
    std::cout << "node_bc_types" << std::endl;
    node_bc_type.print();
    std::cout << "node_bc_values" << std::endl;
    node_bc_vals.print();
  }

  // Get x and y position of nodes for each element
  for (int i = 0; i < NUM_ELEM; i++) {
    for (int j = 0; j < 2; j++) {
      x_elem(j, i) = node_pos(0, elem_list(j, i));
      y_elem(j, i) = node_pos(1, elem_list(j, i));
    }
  }
  
  // Get element lengths
  for (int i = 0; i < NUM_ELEM; i++) {
    elem_lengths(i) = sqrt(pow(x_elem(1, i) - x_elem(0, i), 2) + pow(y_elem(1, i) - y_elem(0, i), 2));
  }
  
  // Calculate the node angle (current implementation only works for convex shapes (0 <= alpha <= pi)
  for (int i = 0; i < NUM_NODE; i++) {
    int prev_node = (i + NUM_NODE - 1) % NUM_NODE;
    int next_node = (i + NUM_NODE + 1) % NUM_NODE;

    double x_prev = node_pos(0, prev_node);
    double y_prev = node_pos(1, prev_node);
    double x = node_pos(0, i);
    double y = node_pos(1, i);
    double x_next = node_pos(0, next_node);
    double y_next = node_pos(1, next_node);

    double a = sqrt(pow(x - x_prev, 2) + pow(y - y_prev, 2));
    double b = sqrt(pow(x_next - x, 2) + pow(y_next - y, 2));
    double c = sqrt(pow(x_next - x_prev, 2) + pow(y_next - y_prev, 2));

    node_angles(i) = acos((pow(a, 2) + pow(b, 2) - pow(c, 2)) / (2 * a * b));
  }

  // Calculate element tangents and normals
  for (int i = 0; i < NUM_ELEM; i++) {
    elem_tangent(0, i) = (x_elem(1, i) - x_elem(0, i)) / elem_lengths(i);
    elem_tangent(1, i) = (y_elem(1, i) - y_elem(0, i)) / elem_lengths(i);

    elem_normal(0, i) = (y_elem(1, i) - y_elem(0, i)) / elem_lengths(i);
    elem_normal(1, i) = -(x_elem(1, i) - x_elem(0, i)) / elem_lengths(i);
  }

  for (int i = 0; i < NUM_NODE; i++) {
    double x = node_pos(0, i);
    double y = node_pos(1, i);

    for (int l = 0; l < NUM_ELEM; l++) {
      int il_1 = elem_list(0, l);
      int il_2 = elem_list(1, l);

      // Handle analytical integrals later
      if (i == il_1 || i == il_2) {
        continue;
      }

      // Gauss quadrature
      for (int k = 0; k < NUM_QUAD; k++) {
        double zeta = quad_points(k);
        double w = quad_weights(k);

        // Get basis fcns
        std::vector<double> phi({(1 - zeta) / 2.0, (1 + zeta) / 2.0});

        // Get sample point from basis functions
        double xk = x_elem(0, l) * phi[0] + x_elem(1, l) * phi[1];
        double yk = y_elem(0, l) * phi[0] + y_elem(1, l) * phi[1];

        // Calculate functions using points
        double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
        double drdn = (elem_normal(0, l) * (xk - x) - elem_normal(1, l) * (yk - y)) / r;

        double G = -log(r);
        double dGdn = -1 / r * drdn;

        A(i, il_1) = A(i, il_1) + phi[0] * dGdn * elem_lengths(l) * w / 2.0;
        A(i, il_2) = A(i, il_2) + phi[1] * dGdn * elem_lengths(l) * w / 2.0;

        B(i, 2*il_1) = B(i, 2*il_1) + phi[0] * G * elem_lengths(l) * w / 2.0;
        B(i, 2*il_2 + 1) = B(i, 2*il_2 + 1) + phi[1] * G * elem_lengths(l) * w / 2.0;

      }
    }
    A(i, i) = A(i, i) + node_angles(i);
  }

  // Analytical integrals
  for (int l = 0; l < NUM_ELEM; l++) {
    int il_1 = elem_list(0, l);
    int il_2 = elem_list(1, l);

    B(il_1, 2 * il_1) = B(il_1, 2 * il_1) + elem_lengths(l) / 2.0 * (1.5 - log(elem_lengths(l)));
    B(il_2, 2 * il_2 + 1) = B(il_2, 2 * il_2 + 1) + elem_lengths(l) / 2.0 * (1.5 - log(elem_lengths(l)));

    B(il_1, 2 * il_2 + 1) = B(il_1, 2 * il_2 + 1) + elem_lengths(l) / 2.0 * (0.5 - log(elem_lengths(l)));
    B(il_2, 2 * il_1) = B(il_2, 2 * il_1) + elem_lengths(l) / 2.0 * (0.5 - log(elem_lengths(l)));
  }

  // Populate F
  for (int i = 0; i < NUM_NODE; i++) {
    for (int s = 0; s < NUM_SRC; s++) {

      double x = node_pos(0, i);
      double y = node_pos(1, i);

      double xs = source_pos(0, s);
      double ys = source_pos(1, s);

      double r = sqrt(pow(xs - x, 2) + pow(ys - y, 2));
      double G = -log(r);

      RHS(i) = RHS(i) + source_values(s) * G;

    }
  }

  A.write("output/A_corner.dat");
  B.write("output/B_corner.dat");
  RHS.write("output/F_corner.dat");

  Array<double> sladek_bcs(NUM_NODE);
  sladek_bcs.fill(NAN);

  for (int j = 0; j < NUM_NODE; j++) {
    // Sladek and Sladek
    if (node_bc_type(j) == 11) {
      int elem_1 = -1;
      int elem_2 = -1;

      for (int l = 0; l < NUM_ELEM; l++) {
        if (elem_list(0, l) == j) {
          elem_2 = l;
        }
        else if (elem_list(1, l) == j) {
          elem_1 = l;
        }
      }
      double dudt_1;
      double dudt_2;
      if (node_bc_type(elem_list(1, elem_1)) == 12) {
        dudt_1 = (node_bc_vals(0, elem_list(1, elem_1)) - node_bc_vals(1, elem_list(0, elem_1))) / elem_lengths(elem_1);
      }
      else {
        dudt_1 = (node_bc_vals(0, elem_list(1, elem_1)) - node_bc_vals(0, elem_list(0, elem_1))) / elem_lengths(elem_1);
      }

      if (node_bc_type(elem_list(1, elem_2)) == 21) {
        dudt_2 = (node_bc_vals(1, elem_list(1, elem_2)) - node_bc_vals(0, elem_list(0, elem_2))) / elem_lengths(elem_2);
      }
      else {
        dudt_2 = (node_bc_vals(0, elem_list(1, elem_2)) - node_bc_vals(0, elem_list(0, elem_2))) / elem_lengths(elem_2);
      }

      double dudt_1_coeff = (elem_normal(0, elem_1) * elem_normal(0, elem_2) + elem_normal(1, elem_1) * elem_normal(0, elem_2)) / (elem_tangent(0, elem_1) * elem_tangent(0, elem_2) + elem_tangent(1, elem_1) * elem_tangent(1, elem_2));
      double dudt_2_coeff = elem_normal(0, elem_1) * elem_tangent(0, elem_2) + elem_normal(1, elem_1) * elem_tangent(1, elem_2) - (elem_tangent(0, elem_1) * elem_tangent(0, elem_2) + elem_tangent(1, elem_1) * elem_tangent(1, elem_2)) * (elem_normal(0, elem_1) * elem_normal(0, elem_2) + elem_normal(1, elem_1) * elem_normal(1, elem_2)) / (elem_tangent(0, elem_1) * elem_normal(0, elem_2) + elem_tangent(1, elem_1) * elem_normal(1, elem_2));

      double dudn_1 = dudt_1_coeff * dudt_1 + dudt_2_coeff * dudt_2;

      std::cout << j << ": " << elem_1 << ", " << elem_2 << ", " << dudt_1 << ", " << dudt_2 << ", " << dudt_1_coeff << ", " << dudt_2_coeff << ", " << dudn_1 << std::endl;
      sladek_bcs(j) = dudn_1;

      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = B(i, 2 * j + 1);
        RHS(i) = RHS(i) - A(i, j) * node_bc_vals(0, j) + B(i, 2 * j) * dudn_1;
      }
    }
    else if (node_bc_type(j) == 12) {
      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = B(i, 2 * j);
        RHS(i) = RHS(i) - A(i, j) * node_bc_vals(0, j) + B(i, 2 * j + 1) * node_bc_vals(1, j);
      }
    }
    else if (node_bc_type(j) == 21) {
      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = B(i, 2 * j + 1);
        RHS(i) = RHS(i) - A(i, j) * node_bc_vals(1, j) + B(i, 2 * j) * node_bc_vals(0, j);
      }

    }
    else if (node_bc_type(j) == 22) {
      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = A(i, j);
        RHS(i) = RHS(i) + B(i, 2 * j) * node_bc_vals(0, j) + B(i, 2 * j + 1) * node_bc_vals(1, j);
      }
    }
    else {
      throw std::invalid_argument("BC type given for node " + std::to_string(j) + " is invalid: " + std::to_string(node_bc_type(j)));
    }
  }

  sladek_bcs.print();

  LHS.write("output/LHS_corner.dat");
  RHS.write("output/RHS_corner.dat");

  lapack_int info;
  lapack_int *ipiv = new lapack_int[LHS.get_rows()];
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, LHS.get_rows(), 1, LHS.dataPtr(),
                       LHS.get_rows(), ipiv, RHS.dataPtr(), RHS.get_rows());
  
  Array<double> u_boundary(NUM_NODE);
  Array<double> dudn_boundary(2, NUM_NODE);

  for (int i = 0; i < NUM_NODE; i++) {
    if (node_bc_type(i) == 11) {
      u_boundary(i) = node_bc_vals(0, i);
      dudn_boundary(0, i) = sladek_bcs(i);
      dudn_boundary(1, i) = RHS(i);
    }
    else if (node_bc_type(i) == 12) {
      u_boundary(i) = node_bc_vals(0, i);
      dudn_boundary(0, i) = RHS(i);
      dudn_boundary(1, i) = node_bc_vals(1, i);
    }
    else if (node_bc_type(i) == 21) {
      u_boundary(i) = node_bc_vals(1, i);
      dudn_boundary(0, i) = node_bc_vals(0, i);
      dudn_boundary(1, i) = RHS(i);
    }
    else if (node_bc_type(i) == 22) {
      u_boundary(i) = RHS(i);
      dudn_boundary(0, i) = node_bc_vals(0, i);
      dudn_boundary(1, i) = node_bc_vals(1, i);
    }
    else {
      std::cout << "bongo";
    }
  }

  u_boundary.write("output/u_boundary_corner.dat");
  dudn_boundary.write("output/dudn_boundary_corner.dat");
}
