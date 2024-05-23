#include "hw3.hpp"
#define NUM_INNER_NODE 24
#define NUM_OUTER_NODE 48
#define NUM_NODE NUM_INNER_NODE + NUM_OUTER_NODE
#define NUM_INNER_ELEM 24
#define NUM_OUTER_ELEM 48
#define NUM_ELEM NUM_INNER_ELEM + NUM_OUTER_ELEM
#define NUM_QUAD 10
#define NUM_SAMP 3816
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433
#define STRUVE_ORDER 20
#define k1 13.84
#define k2 9.28

int main(int argc, char **argv) {
    Array<double> node_pos(2, NUM_NODE);
    Array<int> bc_type(NUM_NODE);
    Array<double> bc_val(NUM_NODE);
    
    Array<int> elem_list(2, NUM_ELEM);
    Array<int> elem_boundary(2, NUM_ELEM);
    
    Array<double> quad_points(NUM_QUAD);
    Array<double> quad_weights(NUM_QUAD);

    Array<double> sample_pos(2, NUM_SAMP);
    Array<int> sample_region(NUM_SAMP);

    Array<double> elem_lengths(NUM_ELEM);

    Array<std::complex<double>> A1(NUM_INNER_NODE, NUM_INNER_NODE);
    Array<std::complex<double>> B1(NUM_INNER_NODE, NUM_INNER_NODE);

    Array<std::complex<double>> A2(NUM_NODE, NUM_NODE);
    Array<std::complex<double>> B2(NUM_NODE, NUM_NODE);

    Array<std::complex<double>> LHS(2 * NUM_INNER_NODE + NUM_OUTER_NODE, 2 * NUM_INNER_NODE + NUM_OUTER_NODE);
    Array<std::complex<double>> LHS_big(4 * NUM_INNER_NODE + NUM_OUTER_NODE, 4 * NUM_INNER_NODE + NUM_OUTER_NODE);
    Array<std::complex<double>> RHS(2 * NUM_INNER_NODE + NUM_OUTER_NODE);
    Array<std::complex<double>> RHS_big(4 * NUM_INNER_NODE + NUM_OUTER_NODE);

    A1.fill(std::complex<double>(0.0, 0.0));
    A2.fill(std::complex<double>(0.0, 0.0));
    B1.fill(std::complex<double>(0.0, 0.0));
    B2.fill(std::complex<double>(0.0, 0.0));
    LHS.fill(std::complex<double>(0.0, 0.0));
    LHS_big.fill(std::complex<double>(0.0, 0.0));
    RHS.fill(std::complex<double>(0.0, 0.0));
    RHS_big.fill(std::complex<double>(0.0, 0.0));

    Array<std::complex<double>> u1(NUM_INNER_NODE);
    Array<std::complex<double>> q1(NUM_INNER_NODE);
    Array<std::complex<double>> u2(NUM_NODE);
    Array<std::complex<double>> q2(NUM_NODE);

    Array<std::complex<double>> u_sample(NUM_SAMP);
    u_sample.fill(std::complex<double>(0.0, 0.0));

    std::string node_filename = "problem_definition/hw3.nod";
    std::string bc_filename = "problem_definition/hw3.bcs";
    std::string elem_filename = "problem_definition/hw3.ele";
    std::string quad_filename =
    "quadrature/deg_" + std::to_string(NUM_QUAD) + ".gqd";
    std::string sample_filename = "problem_definition/hw3.spl";

    std::ifstream node_file(node_filename);
    if (!node_file.is_open()) {
        throw std::invalid_argument("Node file (" + node_filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < NUM_NODE; i++) {
        std::string line;
        std::string token;
        getline(node_file, line);
        std::istringstream line_stream(line);

        // Node number, do nothing
        getline(line_stream, token, ',');

        // Node positions
        getline(line_stream, token, ',');
        node_pos(0, i) = std::stod(token);

        getline(line_stream, token, ',');
        node_pos(1, i) = std::stod(token);
    }

    node_file.close();

    std::ifstream bc_file(bc_filename);
    if (!bc_file.is_open()) {
        throw std::invalid_argument("bc file (" + bc_filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < NUM_NODE; i++) {
        std::string line;
        std::string token;
        getline(bc_file, line);
        std::istringstream line_stream(line);

        // bc number, do nothing
        getline(line_stream, token, ',');

        getline(line_stream, token, ',');
        bc_type(i) = std::stoi(token);

        getline(line_stream, token, ',');
        bc_val(i) = std::stod(token);

    }

    bc_file.close();

    std::ifstream elem_file(elem_filename);
    if (!elem_file.is_open()) {
        throw std::invalid_argument("elem file (" + elem_filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < NUM_ELEM; i++) {
        std::string line;
        std::string token;
        getline(elem_file, line);
        std::istringstream line_stream(line);

        // elem number, do nothing
        getline(line_stream, token, ',');

        getline(line_stream, token, ',');
        elem_list(0, i) = std::stoi(token);

        getline(line_stream, token, ',');
        elem_list(1, i) = std::stoi(token);

        getline(line_stream, token, ',');
        elem_boundary(0, i) = std::stoi(token);

        getline(line_stream, token, ',');
        elem_boundary(1, i) = std::stoi(token);
    }

    elem_file.close();

    std::ifstream quad_file(quad_filename);
    if (!quad_file.is_open()) {
        throw std::invalid_argument("quad file (" + quad_filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < NUM_QUAD; i++) {
        std::string line;
        std::string token;
        getline(quad_file, line);
        std::istringstream line_stream(line);

        // quad number, do nothing
        getline(line_stream, token, ',');

        // Quad point
        getline(line_stream, token, ',');
        quad_points(i) = std::stod(token);

        // Quad weight
        getline(line_stream, token, ',');
        quad_weights(i) = std::stod(token);
    }

    quad_file.close();

    std::ifstream sample_file(sample_filename);
    if (!sample_file.is_open()) {
        throw std::invalid_argument("sample file (" + sample_filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < NUM_SAMP; i++) {
        std::string line;
        std::string token;
        getline(sample_file, line);
        std::istringstream line_stream(line);

        // sample number, do nothing
        getline(line_stream, token, ',');

        // sample positions
        getline(line_stream, token, ',');
        sample_pos(0, i) = std::stod(token);

        getline(line_stream, token, ',');
        sample_pos(1, i) = std::stod(token);

        getline(line_stream, token, ',');
        sample_region(i) = std::stoi(token);
    }

    sample_file.close();

    if (PRINT_INPUTS) {
        std::cout << "Node positions" << std::endl;
        node_pos.print();

        std::cout << "Node BC type" << std::endl;
        bc_type.print();

        std::cout << "Node BC value" << std::endl;
        bc_val.print();

        std::cout << "Element list" << std::endl;
        elem_list.print();

        std::cout << "Element boundary" << std::endl;
        elem_boundary.print();

        std::cout << "Quadrature points" << std::endl;
        quad_points.print();

        std::cout << "Quadrature weights" << std::endl;
        quad_weights.print();

        std::cout << "Sample positions" << std::endl;
        sample_pos.print();

        std::cout << "Sample region" << std::endl;
        sample_region.print();
    }

    for (int i = 0; i < NUM_ELEM; i++) {
        int il1 = elem_list(0, i);
        int il2 = elem_list(1, i);

        double x1 = node_pos(0, il1);
        double x2 = node_pos(0, il2);
        double y1 = node_pos(1, il1);
        double y2 = node_pos(1, il2);

        elem_lengths(i) = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    }

    // Populate A1 and B1 (inner region 1)
    for (int i = 0; i < NUM_INNER_NODE; i++) {
        double x = node_pos(0, i);
        double y = node_pos(1, i);

        double x_behind = node_pos(0, (i + NUM_INNER_NODE - 1) % NUM_INNER_NODE);
        double y_behind = node_pos(1, (i + NUM_INNER_NODE - 1) % NUM_INNER_NODE);
        double x_ahead = node_pos(0, (i + 1) % NUM_INNER_NODE);
        double y_ahead = node_pos(1, (i + 1) % NUM_INNER_NODE);

        double alpha = calc_alpha(x_behind, y_behind, x, y, x_ahead, y_ahead);

        for (int l = 0; l < NUM_INNER_ELEM; l++) {
            int il1, il2;
            if (elem_boundary(0, l) == 1) {
                il1 = elem_list(0, l);
                il2 = elem_list(1, l);
            }
            else {
                throw std::logic_error("Element " + std::to_string(l) + " is either backwards or has the boundary incorrectly described for region 1.\n");
            }

            std::vector<double> x_elem({node_pos(0, il1), node_pos(0, il2)});
            std::vector<double> y_elem({node_pos(1, il1), node_pos(1, il2)});
            std::vector<double> elem_normal({(node_pos(1, il2) - node_pos(1, il1)) / elem_lengths(l), -(node_pos(0, il2) - node_pos(0, il1)) / elem_lengths(l)});

            for (int k = 0; k < NUM_QUAD; k++) {
                double zeta = quad_points(k);
                double w = quad_weights(k);

                if (i == il1 || i == il2) {
                    double a = 0.25;
                    double b = il1 == i ? 0.75 : -0.75;
                    double c = 0.75;
                    double d = -b;

                    double nu = a * pow(zeta, 3) + b * pow(zeta, 2) + c * zeta + d;
                    double dnu = 3 * a * pow(zeta, 2) + 2 * b * zeta + c;

                    std::vector<double> phi({(1 - nu) / 2, (1 + nu) / 2});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    std::complex<double> G = std::complex<double>(0.0, 1.0) / 4.0 * hankel1(0, k1 * r);

                    B1(i, il1) = B1(i, il1) + w * elem_lengths(l) / 2.0 * phi[0] * G * dnu;
                    B1(i, il2) = B1(i, il2) + w * elem_lengths(l) / 2.0 * phi[1] * G * dnu;
                }
                else {
                    std::vector<double> phi({(1 - zeta) / 2.0, (1 + zeta) / 2.0});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    double drdn = (elem_normal[0] * (xk - x) + elem_normal[1] * (yk - y)) / r;
                    std::complex<double> G = std::complex<double>(0, 1.0) / 4.0 * hankel1(0, k1 * r);
                    std::complex<double> dGdr = std::complex<double>(0, -1.0) * k1 / 4.0 * hankel1(1, k1 * r);
                    std::complex<double> dGdn = dGdr * drdn;

                    std::complex<double> a1 = phi[0] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> a2 = phi[1] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b1 = phi[0] * G * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b2 = phi[1] * G * elem_lengths(l) / 2.0 * w;

                    A1(i, il1) = A1(i, il1) + a1;
                    A1(i, il2) = A1(i, il2) + a2;
                    B1(i, il1) = B1(i, il1) + b1;
                    B1(i, il2) = B1(i, il2) + b2;
                }
            }
        }
        A1(i, i) = A1(i, i) + alpha / (2.0 * PI);
    }

    A1.write("output/A1.mat");
    B1.write("output/B1.mat");

    // Populate A2 and B2
    for (int i = 0; i < NUM_NODE; i++) {
        double x = node_pos(0, i);
        double y = node_pos(1, i);
        double alpha;

        // Must handle alpha being switched for boundary between elements
        if (i < NUM_INNER_NODE) {
            double x_behind = node_pos(0, (i + 1) % NUM_INNER_NODE);
            double y_behind = node_pos(1, (i + 1) % NUM_INNER_NODE);
            double x_ahead = node_pos(0, (i + NUM_INNER_NODE - 1) % NUM_INNER_NODE);
            double y_ahead = node_pos(1, (i + NUM_INNER_NODE - 1) % NUM_INNER_NODE);

            alpha = calc_alpha(x_behind, y_behind, x, y, x_ahead, y_ahead);
        }
        else {
            double x_behind = node_pos(0, (i + NUM_OUTER_NODE - NUM_INNER_NODE - 1) % NUM_OUTER_NODE + NUM_INNER_NODE);
            double y_behind = node_pos(1, (i + NUM_OUTER_NODE - NUM_INNER_NODE - 1) % NUM_OUTER_NODE + NUM_INNER_NODE);
            double x_ahead = node_pos(0, (i + 1 - NUM_INNER_NODE) % NUM_OUTER_NODE + NUM_INNER_NODE);
            double y_ahead = node_pos(1, (i + 1 - NUM_INNER_NODE) % NUM_OUTER_NODE + NUM_INNER_NODE);
            alpha = calc_alpha(x_behind, y_behind, x, y, x_ahead, y_ahead);
        }

        for (int l = 0; l < NUM_ELEM; l++) {
            int il1, il2;
            if (elem_boundary(0, l) == 2) {
                il1 = elem_list(0, l);
                il2 = elem_list(1, l);
            }
            // Switch direction of element to have normal be on the right
            else if (elem_boundary(1, l) == 2) {
                il1 = elem_list(1, l);
                il2 = elem_list(0, l);
            }
            else {
                throw std::logic_error("Element " + std::to_string(l) + " has its boundaries defined incorrectly.\n");
            }

            std::vector<double> x_elem({node_pos(0, il1), node_pos(0, il2)});
            std::vector<double> y_elem({node_pos(1, il1), node_pos(1, il2)});
            std::vector<double> elem_normal({(node_pos(1, il2) - node_pos(1, il1)) / elem_lengths(l), -(node_pos(0, il2) - node_pos(0, il1)) / elem_lengths(l)});

            for (int k = 0; k < NUM_QUAD; k++) {
                double zeta = quad_points(k);
                double w = quad_weights(k);

                if (i == il1 || i == il2) {
                    double a = 0.25;
                    double b = il1 == i ? 0.75 : -0.75;
                    double c = 0.75;
                    double d = -b;

                    double nu = a * pow(zeta, 3) + b * pow(zeta, 2) + c * zeta + d;
                    double dnu = 3 * a * pow(zeta, 2) + 2 * b * zeta + c;

                    std::vector<double> phi({(1 - nu) / 2, (1 + nu) / 2});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    std::complex<double> G = std::complex<double>(0.0, 1.0) / 4.0 * hankel1(0, k2 * r);

                    B2(i, il1) = B2(i, il1) + w * elem_lengths(l) / 2.0 * phi[0] * G * dnu;
                    B2(i, il2) = B2(i, il2) + w * elem_lengths(l) / 2.0 * phi[1] * G * dnu;
                }
                else {
                    std::vector<double> phi({(1 - zeta) / 2.0, (1 + zeta) / 2.0});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    double drdn = (elem_normal[0] * (xk - x) + elem_normal[1] * (yk - y)) / r;
                    std::complex<double> G = std::complex<double>(0, 1.0) / 4.0 * hankel1(0, k2 * r);
                    std::complex<double> dGdr = std::complex<double>(0, -1.0) * k2 / 4.0 * hankel1(1, k2 * r);
                    std::complex<double> dGdn = dGdr * drdn;

                    std::complex<double> a1 = phi[0] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> a2 = phi[1] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b1 = phi[0] * G * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b2 = phi[1] * G * elem_lengths(l) / 2.0 * w;

                    A2(i, il1) = A2(i, il1) + a1;
                    A2(i, il2) = A2(i, il2) + a2;
                    B2(i, il1) = B2(i, il1) + b1;
                    B2(i, il2) = B2(i, il2) + b2;
                }
            }
        }
        A2(i, i) = A2(i, i) + alpha / (2.0 * PI);
    }

    A2.write("output/A2.mat");
    B2.write("output/B2.mat");

    // Populate LHS
    // aa portion s
    for (int i = 0; i < NUM_INNER_NODE; i++) {
        for (int j = 0; j < NUM_INNER_NODE; j++) {
            LHS(i, j) = A1(i, j);
            LHS(i, j + NUM_INNER_NODE) = -B1(i, j);
            LHS(i + NUM_INNER_NODE, j) = A2(i, j);
            LHS(i + NUM_INNER_NODE, j + NUM_INNER_NODE) = B2(i, j);
        }
    }

    // ab portion
    for (int i = 0; i < NUM_INNER_NODE; i++) {
        for (int j = 0; j < NUM_OUTER_NODE; j++) {
            LHS(NUM_INNER_NODE + i, 2 * NUM_INNER_NODE + j) = -B2(i, NUM_INNER_NODE + j);
        }
    }

    // ba portions
    for (int i = 0; i < NUM_OUTER_NODE; i++) {
        for (int j = 0; j < NUM_INNER_NODE; j++) {
            LHS(2 * NUM_INNER_NODE + i, j) = A2(NUM_INNER_NODE + i, j);
            LHS(2 * NUM_INNER_NODE + i, NUM_INNER_NODE + j) = B2(NUM_INNER_NODE + i, j);
        }
    }

    // bb portion
    for (int i = 0; i < NUM_OUTER_NODE; i++) {
        for (int j = 0; j < NUM_OUTER_NODE; j++) {
            LHS(2 * NUM_INNER_NODE + i, 2 * NUM_INNER_NODE + j) = - B2(NUM_INNER_NODE + i, NUM_INNER_NODE + j);
        }
    }

    // Populate RHS
    // ab portion
    for (int i = 0; i < NUM_INNER_NODE; i++) {
        for (int j = 0; j < NUM_OUTER_NODE; j++) {
            RHS(NUM_INNER_NODE + i) = RHS(NUM_INNER_NODE + i) - A2(i, NUM_INNER_NODE + j) * bc_val(NUM_INNER_NODE + j);
        }
    }

    // bb portion
    for (int i = 0; i < NUM_OUTER_NODE; i++) {
        for (int j = 0; j < NUM_OUTER_NODE; j++) {
            RHS(2 * NUM_INNER_NODE + i) = RHS(2 * NUM_INNER_NODE + i) - A2(NUM_INNER_NODE + i, NUM_INNER_NODE + j) * bc_val(NUM_INNER_NODE + j);
        }
    }

    LHS.write("output/LHS.mat");
    RHS.write("output/RHS.mat");

    lapack_int info;
    lapack_int *ipiv = new lapack_int[LHS.get_rows()];
    info = LAPACKE_zgesv(LAPACK_COL_MAJOR, LHS.get_rows(), 1, reinterpret_cast <__complex__ double*> (LHS.dataPtr()), LHS.get_rows(), ipiv, reinterpret_cast <__complex__ double*> (RHS.dataPtr()), RHS.get_rows());

    if (info != 0) {
        throw std::logic_error("Matrix failed to solve with code: " + std::to_string(info));
    }

    for (int i = 0; i < NUM_INNER_NODE; i++) {
        u1(i) = RHS(i);
        u2(i) = RHS(i);
        q1(i) = RHS(NUM_INNER_NODE + i);
        q2(i) = -RHS(NUM_INNER_NODE + i);
    }

    for (int i = 0; i < NUM_OUTER_NODE; i++) {
        u2(NUM_INNER_NODE + i) = std::complex<double>(bc_val(NUM_INNER_NODE + i), 0.0);
        q2(NUM_INNER_NODE + i) = RHS(2 * NUM_INNER_NODE + i);
    }

    u1.write("output/u1_boundary.dat");
    u2.write("output/u2_boundary.dat");
    q1.write("output/q1_boundary.dat");
    q2.write("output/q2_boundary.dat");

    // Populate sample positions
    for (int s = 0; s < NUM_SAMP; s++) {
        double xs = sample_pos(0, s);
        double ys = sample_pos(1, s);
        int reg = sample_region(s);
        double k_val = reg == 1 ? k1 : k2;

        for (int l = 0; l < NUM_ELEM; l++) {
            double il1, il2;
            // Outward pointing normal
            if (elem_boundary(0, l) == reg) {
                il1 = elem_list(0, l);
                il2 = elem_list(1, l);
            }
            // Inward pointing normal (switch elem direction)
            else if (elem_boundary(1, l) == reg) {
                il1 = elem_list(1, l);
                il2 = elem_list(0, l);
            }
            // Element does not border region
            else {
                continue;
            }
            std::vector<double> x_elem({node_pos(0, il1), node_pos(0, il2)});
            std::vector<double> y_elem({node_pos(1, il1), node_pos(1, il2)});
            std::vector<double> elem_normal({(node_pos(1, il2) - node_pos(1, il1)) / elem_lengths(l), -(node_pos(0, il2) - node_pos(0, il1)) / elem_lengths(l)});

            // Gauss quadrature
            for (int k = 0; k < NUM_QUAD; k++) {
                double zeta = quad_points(k);
                double w = quad_weights(k);

                std::vector<double> phi({(1 - zeta) / 2.0, (1 + zeta) / 2.0});

                double xk = 0.0;
                double yk = 0.0;

                for (int j = 0; j < 2; j++) {
                    xk += phi[j] * x_elem[j];
                    yk += phi[j] * y_elem[j];
                }

                double r = sqrt(pow(xk - xs, 2) + pow(yk - ys, 2));
                double drdn = (elem_normal[0] * (xk - xs) + elem_normal[1] * (yk - ys)) / r;

                std::complex<double> G = std::complex<double>(0.0, 1.0) / 4.0 * hankel1(0, k_val * r);
                std::complex<double> dGdr = std::complex<double>(0.0, -1.0) * k_val / 4.0 * hankel1(1, k_val * r);
                std::complex<double> dGdn = dGdr * drdn;

                std::vector<std::complex<double>> u(2);
                std::vector<std::complex<double>> q(2);

                if (reg == 1) {
                    u[0] = u1(il1);
                    u[1] = u1(il2);
                    q[0] = q1(il1);
                    q[1] = q1(il2);
                }
                else {
                    u[0] = u2(il1);
                    u[1] = u2(il2);
                    q[0] = q2(il1);
                    q[1] = q2(il2);
                }
                u_sample(s) = u_sample(s) - u[0] * phi[0] * dGdn * elem_lengths(l) / 2.0 * w;
                u_sample(s) = u_sample(s) - u[1] * phi[1] * dGdn * elem_lengths(l) / 2.0 * w;

                u_sample(s) = u_sample(s) + q[0] * phi[0] * G * elem_lengths(l) / 2.0 * w;
                u_sample(s) = u_sample(s) + q[1] * phi[1] * G * elem_lengths(l) / 2.0 * w;
            }
        }
    }

    u_sample.write("output/u_sample.dat");

}

double calc_alpha(double x_behind, double y_behind, double x, double y, double x_ahead, double y_ahead) {
  std::vector<double> A({x_behind - x, y_behind - y});
  std::vector<double> B({x_ahead - x, y_ahead - y});

  double theta = acos((A[0] * B[0] + A[1] * B[1]) / sqrt((pow(A[0], 2) + pow(A[1], 2)) * (pow(B[0], 2) + pow(B[1], 2))));

  double cross = A[0] * B[1] - A[1] * B[0];

  return cross <= 0 ? theta : 2.0 * PI - theta;

}

std::complex<double> hankel1(double nu, double x) {
    double J = std::cyl_bessel_j(nu, x);
    double Y = std::cyl_neumann(nu, x);

    return std::complex<double>(J, Y);
}
