#include "hw1.hpp"
#define NUM_NODE 24
#define NUM_ELEM 24
#define NUM_QUAD 5
#define NUM_SRC 1
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433

int main(int argc, char **argv) {
  // Matrix and vector definitions
  Array2D<double> node_pos(2, NUM_NODE);

  Array2D<int> elem_list(2, NUM_ELEM);

  std::vector<int> bc_types(NUM_NODE);
  std::vector<double> bc_values(NUM_NODE);

  Array2D<double> src_pos(2, NUM_SRC);
  std::vector<double> src_values(NUM_SRC);

  Array2D<double> A(NUM_NODE, NUM_NODE);
  A.fill(0.0);
  Array2D<double> B(NUM_NODE, NUM_NODE);
  B.fill(0.0);

  Array2D<double> LHS(NUM_NODE, NUM_NODE);
  LHS.fill(0.0);

  std::vector<double> F(NUM_NODE);
  std::fill(F.begin(), F.end(), 0.0);

  std::vector<double> RHS(NUM_NODE);
  std::fill(RHS.begin(), RHS.end(), 0.0);

  std::vector<double> quad_points(NUM_QUAD);
  std::vector<double> quad_weights(NUM_QUAD);

  Array2D<double> x_elem(2, NUM_ELEM);
  Array2D<double> y_elem(2, NUM_ELEM);

  std::vector<double> delta_s(NUM_ELEM);

  double alpha = 11 * PI / 24.0;

  // Filename definitions
  std::string pos_filename = "problem_definition/hw1.nod";
  std::string elem_filename = "problem_definition/hw1.ele";
  std::string bc_filename = "problem_definition/hw1.bcs";
  std::string src_filename = "problem_definition/hw1.src";
  std::string quad_filename =
      "quadrature/deg_" + std::to_string(NUM_QUAD) + ".gqd";

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
    node_pos(0, i) = std::stod(token);

    // Node y
    getline(lineStream, token, ',');
    node_pos(1, i) = std::stod(token);
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
    elem_list(0, i) = std::stoi(token);

    // Second node
    getline(lineStream, token, ',');
    elem_list(1, i) = std::stoi(token);
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
    src_pos(0, i) = std::stod(token);

    // Source y pos
    getline(lineStream, token, ',');
    src_pos(1, i) = std::stod(token);

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

  for (int i = 0; i < NUM_QUAD; i++) {
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
    node_pos.print();

    std::cout << "Element list: " << std::endl;
    elem_list.print();

    std::cout << "BC types: " << std::endl;
    print_vector(bc_types);

    std::cout << "BC values: " << std::endl;
    print_vector(bc_values);

    std::cout << "Source positions: " << std::endl;
    src_pos.print();

    std::cout << "Source values: " << std::endl;
    print_vector(src_values);

    std::cout << "Quadrature points: " << std::endl;
    print_vector(quad_points);

    std::cout << "Quadrature weights: " << std::endl;
    print_vector(quad_weights);
  }

  // Fill x_elem and y_elem
  for (int l = 0; l < NUM_ELEM; l++) {
    for (int k = 0; k < 2; k++) {
      x_elem(k, l) = node_pos(0, elem_list(k, l));
      y_elem(k, l) = node_pos(1, elem_list(k, l));
    }
  }

  // Fill delta_s
  for (int l = 0; l < NUM_ELEM; l++) {
    delta_s[l] = sqrt(pow(x_elem(1, l) - x_elem(0, l), 2.0) +
                      pow(y_elem(1, l) - y_elem(0, l), 2));
  }

  // Populate A and B
  for (int i = 0; i < NUM_NODE; i++) {
    for (int l = 0; l < NUM_ELEM; l++) {
      if (i == elem_list(0, l) || i == elem_list(1, l)) {
        A(elem_list(0, l), elem_list(0, l)) += alpha;
        A(elem_list(1, l), elem_list(1, l)) += alpha;

        B(elem_list(0, l), elem_list(0, l)) +=
            delta_s[l] * (1.5 - log(delta_s[l])) / 2.0;
        B(elem_list(0, l), elem_list(1, l)) +=
            delta_s[l] * (0.5 - log(delta_s[l])) / 2.0;
        B(elem_list(1, l), elem_list(0, l)) +=
            delta_s[l] * (1.5 - log(delta_s[l])) / 2.0;
        B(elem_list(1, l), elem_list(1, l)) +=
            delta_s[l] * (0.5 - log(delta_s[l])) / 2.0;

        continue;
      }

      for (int k = 0; k < NUM_QUAD; k++) {
        double zeta = quad_points[k];
        double w = quad_weights[k];
        std::vector<double> phi{(1 - zeta) / 2.0, (1 + zeta) / 2.0};
        double xs = 0;
        double ys = 0;

        for (int j = 0; j < 2; j++) {
          xs += x_elem(j, l) * phi[j];
          ys += y_elem(j, l) * phi[j];
        }

        double r =
            sqrt(pow(xs - node_pos(0, i), 2) + pow(ys - node_pos(1, i), 2));

        double drdn = ((y_elem(1, l) - y_elem(0, l)) * (xs - node_pos(0, i)) -
                       (x_elem(1, l) - x_elem(0, l)) * (ys - node_pos(1, i))) /
                      (delta_s[l] * r);

        double G = -log(r);
        double dGdn = -drdn / r;

        for (int j = 0; j < 2; j++) {
          A(i, elem_list(j, l)) += w * phi[j] * dGdn * delta_s[l] / 2.0;
          B(i, elem_list(j, l)) += w * phi[j] * G * delta_s[l] / 2.0;
        }
      }
    }
    A(i, i) += alpha;
  }

  // Populate F
  for (int i = 0; i < NUM_NODE; i++) {
    for (int s = 0; s < NUM_SRC; s++) {
      double r = sqrt(pow(src_pos(0, s) - node_pos(0, i), 2) +
                      pow(src_pos(1, s) - node_pos(1, i), 2));

      F[i] += -src_values[s] * log(r);
    }
  }

  A.write("output/A.dat");
  B.write("output/B.dat");
  write_vector(F, "output/F.dat");

  // Populate LHS and RHS
  for (int i = 0; i < NUM_NODE; i++) {
    if (bc_types[i] == 1) {
      F -= A.getColVec(i) * bc_values[i];
      LHS.insert(-B, 0, i, NUM_NODE, 1, 0, i);
    } else {
      F += B.getColVec(i) * bc_values[i];
      LHS.insert(A, 0, i, NUM_NODE, 1, 0, i);
    }
  }

  LHS.write("output/LHS.dat");
  write_vector(F, "output/RHS.dat");

  lapack_int info;
  lapack_int *ipiv = new lapack_int[A.getRows()];
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, LHS.getRows(), 1, LHS.dataPtr(),
                       A.getRows(), ipiv, F.data(), F.size());

  write_vector(F, "output/new_bcs.dat");
  return EXIT_SUCCESS;
}

template <typename T> void print_vector(std::vector<T> vec) {
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    std::cout << std::setw(10) << std::setprecision(4) << *it << std::endl;
  }
}

template <typename T>
void write_vector(std::vector<T> vec, std::string filename) {
  std::ofstream file;
  file.open(filename);
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    file << *it << std::endl;
  }
  file.close();
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  if (a.size() != b.size()) {
    std::cerr << "Vectors a and b in a+b are not the same size" << std::endl;
    return a;
  } else {
    std::vector<T> vec3(a.size());
    for (int i = 0; i < a.size(); i++) {
      vec3[i] = a[i] + b[i];
    }
    return vec3;
  }
}

template <typename T>
void operator+=(std::vector<T> &a, const std::vector<T> &b) {
  if (a.size() != b.size()) {
    std::cerr << "Vectors a and b in a+b are not the same size" << std::endl;
  } else {
    for (int i = 0; i < a.size(); i++) {
      a[i] += b[i];
    }
  }
}

template <typename T>
void operator-=(std::vector<T> &a, const std::vector<T> &b) {
  if (a.size() != b.size()) {
    std::cerr << "Vectors a and b in a-b are not the same size" << std::endl;
  } else {
    for (int i = 0; i < a.size(); i++) {
      a[i] -= b[i];
    }
  }
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &vec, const T &scalar) {
  std::vector<T> out(vec.size());
  for (std::size_t i = 0; i < vec.size(); i++) {
    out[i] = vec[i] * scalar;
  }
  return out;
}

template <typename T>
std::vector<T> operator*(const T &scalar, const std::vector<T> &vec) {
  std::vector<T> out(vec.size());
  for (std::size_t i = 0; i < vec.size(); i++) {
    out[i] = vec[i] * scalar;
  }
  return out;
}
