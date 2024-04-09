#include "hw1.hpp"
#include <lapacke.h>
#define NUM_NODE 24
#define NUM_ELEM 24
#define NUM_QUAD 5
#define NUM_SRC 1
#define NUM_SAMP 721
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433

double calc_alpha(double x_behind, double y_behind, double x, double y, double x_ahead, double y_ahead);

int main(int argc, char **argv) {
  // Matrix and vector definitions
  Array2D<double> node_pos(2, NUM_NODE);

  Array2D<int> elem_list(2, NUM_ELEM);

  std::vector<int> bc_types(NUM_NODE);
  std::vector<double> bc_values(NUM_NODE);

  Array2D<double> source_pos(2, NUM_SRC);
  std::vector<double> source_values(NUM_SRC);

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

  std::vector<double> alpha(NUM_NODE);

  std::vector<double> u(NUM_ELEM);
  std::vector<double> dudn(NUM_ELEM);

  Array2D<double> samp_points(2, NUM_SAMP);
  std::vector<double> u_sample(NUM_SAMP);
  std::fill(u_sample.begin(), u_sample.end(), 0.0);

  Array2D<double> u_elem(2, NUM_ELEM);
  Array2D<double> dudn_elem(2, NUM_ELEM);

  // Filename definitions
  std::string pos_filename = "problem_definition/hw1.nod";
  std::string elem_filename = "problem_definition/hw1.ele";
  std::string bc_filename = "problem_definition/hw1.bcs";
  std::string source_filename = "problem_definition/hw1.src";
  std::string quad_filename =
      "quadrature/deg_" + std::to_string(NUM_QUAD) + ".gqd";
  std::string samp_filename = "problem_definition/domain_samples.dat";

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
    source_values[i] = std::stod(token);
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

  // Sample points
  std::ifstream samp_file(samp_filename);

  if (!samp_file.is_open()) {
    std::cerr << "Could not open sample file" << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0; i < NUM_SAMP; i++) {
    std::string line;
    std::string token;
    getline(samp_file, line);
    std::istringstream lineStream(line);

    // Point number, do nothing
    getline(lineStream, token, ',');

    // X coordinate of sample point;
    getline(lineStream, token, ',');
    samp_points(0, i) = std::stod(token);

    // Y coordinate of sample point;
    getline(lineStream, token, ',');
    samp_points(1, i) = std::stod(token);
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
    source_pos.print();

    std::cout << "Source values: " << std::endl;
    print_vector(source_values);

    std::cout << "Quadrature points: " << std::endl;
    print_vector(quad_points);

    std::cout << "Quadrature weights: " << std::endl;
    print_vector(quad_weights);

    std::cout << "Sample points" << std::endl;
    samp_points.print();
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
                      pow(y_elem(1, l) - y_elem(0, l), 2.0));
  }

  // Fill alpha
  for (int i = 0; i < NUM_NODE; i++) {
    double x_ahead = node_pos(0, (i + 1) % NUM_NODE);
    double x = node_pos(0, i);
    double y_ahead = node_pos(1, (i + 1) % NUM_NODE);
    double y = node_pos(1, i);
    double x_behind = node_pos(0, (i - 1) % NUM_NODE);
    double y_behind = node_pos(1, (i - 1) % NUM_NODE);

    if (i == 0) {
      x_behind = node_pos(0, NUM_NODE - 1);
      y_behind = node_pos(1, NUM_NODE - 1);
    }

    alpha[i] = calc_alpha(x_behind, y_behind, x, y, x_ahead, y_ahead);

  }

  // Populate A and B
  for (int i = 0; i < NUM_NODE; i++) {
    for (int l = 0; l < NUM_ELEM; l++) {
      if (i == elem_list(0, l) || i == elem_list(1, l)) {

        //A(elem_list(0, l), elem_list(0, l)) = A(elem_list(0, l), elem_list(0, l)) + alpha[elem_list(0, l)] / 2.0;

        //A(elem_list(1, l), elem_list(1, l)) = A(elem_list(1, l), elem_list(1, l)) + alpha[elem_list(1, l)] / 2.0;

        B(elem_list(0, l), elem_list(0, l)) = B(elem_list(0, l), elem_list(0, l)) + delta_s[l] / 2.0 * (1.5 - log(delta_s[l]));

        B(elem_list(0, l), elem_list(1, l)) = B(elem_list(0, l), elem_list(1, l)) + delta_s[l] / 2.0 * (0.5 - log(delta_s[l]));

        B(elem_list(1, l), elem_list(1, l)) = B(elem_list(1, l), elem_list(1, l)) + delta_s[l] / 2.0 * (1.5 - log(delta_s[l]));

        B(elem_list(1, l), elem_list(0, l)) = B(elem_list(1, l), elem_list(0, l)) + delta_s[l] / 2.0 * (0.5 - log(delta_s[l]));

      }
      else {
        for (int k = 0; k < NUM_QUAD; k++) {
          double zeta = quad_points[k];
          double w = quad_weights[k];
          std::vector<double> phi({(1 - zeta) / 2.0, (1 + zeta) / 2.0});

          double x = node_pos(0, i);
          double y = node_pos(1, i);

          double xs = x_elem(0, l) * phi[0] + x_elem(1, l) * phi[1];
          double ys = y_elem(0, l) * phi[0] + y_elem(1, l) * phi[1];

          double r = sqrt(pow(xs - x, 2) + pow(ys - y, 2));
          double drdn = ((y_elem(1, l) - y_elem(0, l)) * (xs - x) - (x_elem(1, l) - x_elem(0, l)) * (ys - y)) / (delta_s[l] * r);

          double G = -log(r);
          double dGdn = -1 / r * drdn;

          //std::cout << i << "," << l << "," << k << "," << zeta << "," << w << "," << phi[0] << "," << phi[1] << "," << x << "," << y << "," << x_elem(0, l) << "," << y_elem(0, l) << "," << x_elem(1, l) << "," << y_elem(1, l) << "," << xs << "," << ys << "," << delta_s[l] << "," << r << "," << drdn << "," << G << "," << dGdn << std::endl;

          A(i, elem_list(0, l)) = A(i, elem_list(0, l)) + phi[0] * dGdn * delta_s[l] * w / 2.0;
          A(i, elem_list(1, l)) = A(i, elem_list(1, l)) + phi[1] * dGdn * delta_s[l] * w / 2.0;

          B(i, elem_list(0, l)) = B(i, elem_list(0, l)) + phi[0] * G * delta_s[l] * w / 2.0;
          B(i, elem_list(1, l)) = B(i, elem_list(1, l)) + phi[1] * G * delta_s[l] * w / 2.0;

        }
      }
      //std::cout << i << ", " << l << ", " << elem_list(0, l) << ", " << elem_list(1, l) << ", " <<B(1,0) << std::endl;
    }
    A(i,i) = A(i,i) + alpha[i];
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

      F[i] = F[i] + source_values[s] * G;

    }
  }

  A.write("A.dat");
  B.write("B.dat");
  write_vector(F, "F.dat");


  for (int j = 0; j < NUM_NODE; j++) {
    if (bc_types[j] == 1) {
      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = -B(i, j);
        F[i] = F[i] - A(i, j) * bc_values[j];
      }
    }
    else {
      for (int i = 0; i < NUM_NODE; i++) {
        LHS(i, j) = A(i, j);
        F[i] = F[i] + B(i, j) * bc_values[j];
      }
    }
  }

  LHS.write("LHS.dat");
  write_vector(F, "RHS.dat");

  // Solve system
  lapack_int info;
  lapack_int *ipiv = new lapack_int[A.getRows()];
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, LHS.getRows(), 1, LHS.dataPtr(),
                       A.getRows(), ipiv, F.data(), F.size());

  write_vector(F, "out.dat");

  for (int i = 0; i < NUM_NODE; i++) {
    if (bc_types[i] == 1) {
      u[i] = bc_values[i];
      dudn[i] = F[i];
    }
    else {
      u[i] = F[i];
      dudn[i] = bc_values[i];
    }
  }

  write_vector(u, "u.dat");
  write_vector(dudn, "dudn.dat");

}

double calc_alpha(double x_behind, double y_behind, double x, double y, double x_ahead, double y_ahead) {
  double a = sqrt(pow(x_ahead - x, 2) + pow(y_ahead - y, 2));
  double b = sqrt(pow(x - x_behind, 2) + pow(y - y_behind, 2));
  double c = sqrt(pow(x_ahead - x_behind, 2) + pow(y_ahead - y_behind, 2));

  return acos((pow(a, 2) + pow(b, 2) - pow(c, 2)) / (2 * a * b));
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
