#include "hw3.hpp"
#define NUM_NODE 24
#define NUM_ELEM 24
#define NUM_QUAD 64
#define NUM_SAMP 1920
#define PRINT_INPUTS false
#define PI 3.1415926535897932384626433
#define STRUVE_ORDER 20
#define kappa 1.0

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

    Array<std::complex<double>> A(NUM_NODE, NUM_NODE);
    Array<std::complex<double>> B(NUM_NODE, NUM_NODE);

    Array<std::complex<double>> LHS(NUM_NODE, NUM_NODE);
    Array<std::complex<double>> RHS(NUM_NODE);

    A.fill(std::complex<double>(0.0, 0.0));
    B.fill(std::complex<double>(0.0, 0.0));
    LHS.fill(std::complex<double>(0.0, 0.0));
    RHS.fill(std::complex<double>(0.0, 0.0));

    Array<std::complex<double>> u(NUM_NODE);
    Array<std::complex<double>> q(NUM_NODE);
    Array<std::complex<double>> u_sample(NUM_SAMP);
    u_sample.fill(std::complex<double>(0.0, 0.0));

    std::string node_filename = "problem_definition/simple.nod";
    std::string bc_filename = "problem_definition/simple.bcs";
    std::string elem_filename = "problem_definition/simple.ele";
    std::string quad_filename =
    "quadrature/deg_" + std::to_string(NUM_QUAD) + ".gqd";
    std::string sample_filename = "problem_definition/simple.spl";


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

    for (int i = 0; i < NUM_NODE; i++) {
        double x = node_pos(0, i);
        double y = node_pos(1, i);
        double x_behind, y_behind, x_ahead, y_ahead;

        if (i == 0) {
            x_behind = node_pos(0, NUM_NODE - 1);
            y_behind = node_pos(1, NUM_NODE - 1);
            x_ahead = node_pos(0, 1);
            y_ahead = node_pos(1, 1);
        }
        else if (i == NUM_NODE - 1) {
            x_behind = node_pos(0, NUM_NODE - 2);
            y_behind = node_pos(1, NUM_NODE - 2);
            x_ahead = node_pos(0, 0);
            y_ahead = node_pos(1, 0);
        }
        else {
            x_behind = node_pos(0, i - 1);
            y_behind = node_pos(1, i - 1);
            x_ahead = node_pos(0, i + 1);
            y_ahead = node_pos(1, i + 1);
        }

        double alpha = calc_alpha(x_behind, y_behind, x, y, x_ahead, y_ahead);

        for (int l = 0; l < NUM_ELEM; l++) {
            int il1, il2;

            if (elem_boundary(0, l) == 1) {
                il1 = elem_list(0, l);
                il2 = elem_list(1, l);
            }
            else if (elem_boundary(1, l) == 1) {
                il1 = elem_list(1, l);
                il2 = elem_list(0, l);
            }
            else {
                continue;
            }

            std::vector<double> x_elem({node_pos(0, il1), node_pos(0, il2)});
            std::vector<double> y_elem({node_pos(1, il1), node_pos(1, il2)});
            std::vector<double> elem_normal({(y_elem[1] - y_elem[0]) / elem_lengths(l), -(x_elem[1] - x_elem[0]) / elem_lengths(l)});

            for (int k = 0; k < NUM_QUAD; k++) {
                double z = quad_points(k);
                double w = quad_weights(k);

                if (il1 == i || il2 == i) {
                    double a = 0.25;
                    double b = il1 == i ? 0.75 : -0.75;
                    double c = 0.75;
                    double d = -b;

                    double nu = a * pow(z, 3) + b * pow(z, 2) + c * z + d;
                    double dnu = 3 * a * pow(z, 2) + 2 * b * z + c;

                    std::vector<double> phi({(1 - nu) / 2, (1 + nu) / 2});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    std::complex<double> G = std::complex<double>(0.0, 1.0) / 4.0 * hankel1(0, kappa * r);

                    B(i, il1) = B(i, il1) + w * elem_lengths(l) / 2.0 * phi[0] * G * dnu;
                    B(i, il2) = B(i, il2) + w * elem_lengths(l) / 2.0 * phi[1] * G * dnu;
                }
                else {
                    std::vector<double> phi({(1 - z) / 2, (1 + z) / 2});

                    double xk = 0.0;
                    double yk = 0.0;

                    for (int j = 0; j < 2; j++) {
                        xk += phi[j] * x_elem[j];
                        yk += phi[j] * y_elem[j];
                    }

                    double r = sqrt(pow(xk - x, 2) + pow(yk - y, 2));
                    double drdn = (elem_normal[0] * (xk - x) + elem_normal[1] * (yk - y)) / r;
                    std::complex<double> G = std::complex<double>(0, 1.0) / 4.0 * hankel1(0, kappa * r);
                    std::complex<double> dGdr = std::complex<double>(0, -1.0) * kappa / 4.0 * hankel1(1, kappa * r);
                    std::complex<double> dGdn = dGdr * drdn;

                    std::complex<double> a1 = phi[0] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> a2 = phi[1] * dGdn * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b1 = phi[0] * G * elem_lengths(l) / 2.0 * w;
                    std::complex<double> b2 = phi[1] * G * elem_lengths(l) / 2.0 * w;

                    A(i, il1) = A(i, il1) + a1;
                    A(i, il2) = A(i, il2) + a2;
                    B(i, il1) = B(i, il1) + b1;
                    B(i, il2) = B(i, il2) + b2;

                }

            }
        }
        A(i, i) = A(i, i) + alpha / (2.0 * PI);
    }

    A.write("output/A_simple.dat");
    B.write("output/B_simple.dat");

    for (int i = 0; i < NUM_NODE; i++) {
        for (int j = 0; j < NUM_NODE; j++) {
            LHS(i, j) = -B(i, j);
            RHS(i) = RHS(i) - bc_val(j) * A(i, j);
        }
    }

    LHS.write("output/LHS_simple.dat");
    RHS.write("output/RHS_simple.dat");

    lapack_int info;
    lapack_int *ipiv = new lapack_int[LHS.get_rows()];
    info = LAPACKE_zgesv(LAPACK_COL_MAJOR, LHS.get_rows(), 1, reinterpret_cast <__complex__ double*> (LHS.dataPtr()), LHS.get_rows(), ipiv, reinterpret_cast <__complex__ double*> (RHS.dataPtr()), RHS.get_rows());

    if (info != 0) {
        throw std::logic_error("Matrix failed to solve with code: " + std::to_string(info));
    }

    for (int i = 0; i < NUM_NODE; i++) {
        u(i) = bc_val(i);
        q(i) = RHS(i);
    }

    u.write("output/u_simple_boundary.dat");
    q.write("output/q_simple_boundary.dat");

    for (int s = 0; s < NUM_SAMP; s++) {
        double xs = sample_pos(0, s);
        double ys = sample_pos(1, s);

        for (int l = 0; l < NUM_ELEM; l++) {
            int il1, il2;

            if (elem_boundary(0, l) == 1) {
                il1 = elem_list(0, l);
                il2 = elem_list(1, l);
            }
            else if (elem_boundary(1, l) == 1) {
                il1 = elem_list(1, l);
                il2 = elem_list(0, l);
            }
            else {
                continue;
            }

            std::vector<double> x_elem({node_pos(0, il1), node_pos(0, il2)});
            std::vector<double> y_elem({node_pos(1, il1), node_pos(1, il2)});
            std::vector<double> elem_normal({(y_elem[1] - y_elem[0]) / elem_lengths(l), -(x_elem[1] - x_elem[0]) / elem_lengths(l)});

            for (int k = 0; k < NUM_QUAD; k++) {
                double z = quad_points(k);
                double w = quad_weights(k);

                std::vector<double> phi({(1 - z) / 2, (1 + z) / 2});

                double xk = 0.0;
                double yk = 0.0;

                for (int j = 0; j < 2; j++) {
                    xk += phi[j] * x_elem[j];
                    yk += phi[j] * y_elem[j];
                }

                double r = sqrt(pow(xk - xs, 2) + pow(yk - ys, 2));
                double drdn = (elem_normal[0] * (xk - xs) + elem_normal[1] * (yk - ys)) / r;

                std::complex<double> G = std::complex<double>(0.0, 1.0) / 4.0 * hankel1(0, kappa * r);
                std::complex<double> dGdr = std::complex<double>(0.0, -1.0) * kappa / 4.0 * hankel1(1, kappa * r);
                std::complex<double> dGdn = dGdr * drdn;

                u_sample(s) = u_sample(s) - u(il1) * phi[0] * dGdn * elem_lengths(l) / 2.0 * w;
                u_sample(s) = u_sample(s) - u(il2) * phi[1] * dGdn * elem_lengths(l) / 2.0 * w;

                u_sample(s) = u_sample(s) + q(il1) * phi[0] * G * elem_lengths(l) / 2.0 * w;
                u_sample(s) = u_sample(s) + q(il2) * phi[1] * G * elem_lengths(l) / 2.0 * w;
            }
        }
    }

    u_sample.write("output/u_simple_sample.dat");

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

int factorial(int n) {
    int out = 1;
    for (int i = 2; i <= n; i++) {
        out *= i;
    }

    return out;
}

int factorial2(int n) {
    int out = 1;
    for (int i = n; i > 1; i -= 2) {
        out *= i;
    }

    return out;
}

double struve(double nu, double x) {
    double out = 0;
    if (nu == 0) {
        for (int i = 0; i <= STRUVE_ORDER; i++) {
            double val = 2.0 / PI * pow(x, 2*i+1) / pow((double) factorial2(2*i+1), 2);

            out += i % 2 == 0 ? val : -val;
        }
    }
    else if (nu == 1) {
        for (int i = 1; i <= STRUVE_ORDER; i++) {
            double val = 2.0 / PI * pow(x, 2*i) / (pow((double) factorial2(2*i-1), 2) * ((double) 2 * i + 1));

            out += i % 2 == 0 ? -val : val;
        }
    }
    else {
        throw std::logic_error("nu < 1 not implemented yet\n");
    }

    return out;

}

std::complex<double> int_hankel1(double x) {
    return x * hankel1(0, x) + PI / 2.0 * x * (struve(0, x) * hankel1(1, x) - struve(1, x) * hankel1(0, x));
}
