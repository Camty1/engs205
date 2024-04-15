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
  Array<int> node_bc_types(NUM_NODE);
  Array<double> node_bc_coeffs(3, NUM_NODE);
  
  Array<int> elem_list(2, NUM_ELEM);
  Array<int> elem_boundaries(2, NUM_ELEM);
  Array<int> elem_bc_types(NUM_ELEM);
  Array<double> elem_bc_values(NUM_ELEM);

  Array<double> source_pos(2, NUM_SRC);
  Array<double> source_values(NUM_SRC);

  Array<double> quad_points(NUM_QUAD);
  Array<double> quad_weights(NUM_QUAD);

  Array<double> elem_node_bc_flags(2, NUM_ELEM);

  Array<double> x_elem(2, NUM_ELEM);
  Array<double> y_elem(2, NUM_ELEM);

  Array<double> LHS(NUM_NODE, NUM_NODE);
  LHS.fill(0.0);

  Array<double> RHS(NUM_NODE);
  RHS.fill(0.0);

  // Filenames
  std::string node_filename = "problem_definition/hw1.nodb";
  std::string elem_filename = "problem_definition/hw1.eleb";
  std::string source_filename = "problem_definition/hw1.src";
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
    node_bc_types(i) = std::stoi(token);

    // Node BC coeffs
    for (int j = 0; j < 3; j++) {
      getline(lineStream, token, ',');
      node_bc_coeffs(j, i) = std::stod(token);
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

    // Element boundaries
    for (int j = 0; j < 2; j++) {
      getline(lineStream, token, ',');
      elem_boundaries(j, i) = std::stoi(token);
    }

    // 0 is external, which should be to the right of the element
    if (elem_boundaries(0, i) == 0) {
      int temp = elem_list(0, i);
      elem_list(0, i) = elem_list(1, i);
      elem_list(1, i) = temp;

      temp = elem_boundaries(0, i);
      elem_boundaries(0, i) = elem_boundaries(1, i);
      elem_boundaries(1, i) = temp;
    }

    // BC type
    getline(lineStream, token, ',');
    elem_bc_types(i) = std::stoi(token);

    // BC value
    getline(lineStream, token, ',');
    elem_bc_values(i) = std::stod(token);

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
    node_bc_types.print();
    std::cout << "elem_bc_types" << std::endl;
    elem_bc_types.print();
  }

  // Set node element bc flags 1 is type 1, 2 is type 2, 0 is unknown
  for (int i = 0; i < NUM_ELEM; i++) {
    if (node_bc_types(elem_list(0, i)) == 1 && node_bc_types(elem_list(1, i)) == 1) {
      elem_node_bc_flags(0, i) = 01;
      elem_node_bc_flags(1, i) = 01;
    }
    else if (elem_bc_types(i) == 2) {
      for (int j = 0; j < 2; j++) {
        if (node_bc_types(elem_list(j, i)) == 1) {
          elem_node_bc_flags(j, i) = 21;
        }
        else {
          node_bc_types(elem_list(j, i)) = 2;
          node_bc_coeffs(0, elem_list(j, i)) = 0.0;
          node_bc_coeffs(1, elem_list(j, i)) = 1.0;

          elem_node_bc_flags(j, i) = 20;
        }
      }
    }
    else {

      throw std::invalid_argument("BCs are inconsistent on elem " + std::to_string(i) + ", which contains nodes " + std::to_string(elem_list(0, i)) + " and " + std::to_string(elem_list(1, i)) + ".  The BC type of the element is " + std::to_string(elem_bc_types(i)) + " and the node BC types are " + std::to_string(node_bc_types(elem_list(0, i))) + " and " + std::to_string(node_bc_types(elem_list(1, i))) + ".");
    }
  }

  // Get x and y position of nodes for each element
  for (int i = 0; i < NUM_ELEM; i++) {
    for (int j = 0; j < 2; j++) {
      x_elem(j, i) = node_pos(0, elem_list(j, i));
      y_elem(j, i) = node_pos(1, elem_list(j, i));
    }
  }



}
