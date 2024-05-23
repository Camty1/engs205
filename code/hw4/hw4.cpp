#include "hw4.hpp"
#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    int n = 4;
    if (argc == 2) {
        n = std::stod(argv[1]);
    }

    Array<double> x_values(n);
    Array<double> func_values(n);
    Array<double> analytical_deriv_values(n);
    Array<double> finite_diff_values(n);
    Array<int> elem_list(2, n - 1);
    Array<double> elem_mat(2, 2);
    Array<double> elem_vec(2);
    Array<double> fem_values(n + 1);
    Array<std::complex<double>> spec_values(n);
    Array<double> hanning_values(n);
    Array<double> hanning_window_values(n);

    fftw_plan p_forward = fftw_plan_dft_r2c_1d(n, func_values.dataPtr(), (fftw_complex*) spec_values.dataPtr(), FFTW_ESTIMATE);
    fftw_plan p_inverse = fftw_plan_dft_c2r_1d(n, (fftw_complex*) spec_values.dataPtr(), func_values.dataPtr(), FFTW_ESTIMATE);
    fftw_plan p_forward_hanning = fftw_plan_dft_r2c_1d(n, hanning_values.dataPtr(), (fftw_complex*) spec_values.dataPtr(), FFTW_ESTIMATE);
    fftw_plan p_inverse_hanning = fftw_plan_dft_c2r_1d(n, (fftw_complex*) spec_values.dataPtr(), hanning_values.dataPtr(), FFTW_ESTIMATE);

    double dx = 2.0 * PI / (double) (n - 1);

    hanning_window(n, &hanning_window_values);

    for (int i = 0; i < n; i++) {
        x_values(i) = (double) i * dx;
        func_values(i) = f(x_values(i));
        analytical_deriv_values(i) = dfdx(x_values(i));
        spec_values(i) = func_values(i);
        hanning_values(i) = func_values(i);
    }

    analytical_deriv_values.write("output/analytical_" + std::to_string(n) + ".dat");

    for (int i = 1; i < n - 1; i++) {
        finite_diff_values(i) = (func_values(i+1) - func_values(i-1)) / (2.0 * dx);
    }
    finite_diff_values(0) = (func_values(1) - func_values(n-2)) / (2.0 * dx);
    finite_diff_values(n-1) = finite_diff_values(0);

    finite_diff_values.write("output/fd_" + std::to_string(n) + ".dat");

    for (int i = 0; i < n - 1; i++) {
        elem_list(0, i) = i;
        elem_list(1, i) = i + 1;
    }
    int kl = 1;
    int ku = 1; 

    Array<double> A_fem(2*kl + ku + 1, n);
    Array<double> b_fem(n);
    A_fem.fill(0.0);
    b_fem.fill(0.0);

    for (int l = 0; l < n - 1; l++) {
        fem_elem_matrix(dx, func_values(elem_list(0, l)), func_values(elem_list(1, l)), &elem_mat, &elem_vec);

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                int row = elem_list(i, l);
                int col = elem_list(j, l);

                row = kl + ku + row - col;

                A_fem(row, col) = A_fem(row, col) + elem_mat(i, j);

            }
        }

        for (int i = 0; i < 2; i++) {
            int row = elem_list(i, l);
            b_fem(row) = b_fem(row) + elem_vec(i);
        }
    }
    lapack_int info;
    lapack_int *ipiv = new lapack_int[b_fem.get_rows()];
    info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, n, n, kl, ku, A_fem.dataPtr(), 2 * kl + ku + 1, ipiv);

    if (info != 0) {
        throw(std::logic_error("FEM matrix failed to factor"));
    }

    info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', n, kl, ku, 1, A_fem.dataPtr(), 2 * kl + ku + 1, ipiv, b_fem.dataPtr(), n);

    if (info != 0) {
        throw(std::logic_error("FEM matrix failed to solve"));
    }
    b_fem.write("output/fem_" + std::to_string(n) + ".dat");
    
    fftw_execute(p_forward);

    for (int i = 0; i < n; i++) {
        spec_values(i) = spec_values(i) / (double) n;
    }
    for (int k = 0; k < n; k++) {
        if (k < n / 2) {
            spec_values(k) = (std::complex<double>) 1j * (double) k * spec_values(k);
        }
        else {
            spec_values(k) = (std::complex<double>) 1j * (double) (k - n) * spec_values(k);
        }
    }
    fftw_execute(p_inverse);
    func_values.write("output/sm_" + std::to_string(n) + ".dat");

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_inverse);

    fftw_execute(p_forward_hanning);
    for (int i = 0; i < n; i++) {
        spec_values(i) = spec_values(i) / (double) n * hanning_window_values(i);
    }
    for (int k = 0; k < n; k++) {
        if (k < n / 2) {
            spec_values(k) = (std::complex<double>) 1j * (double) k * spec_values(k);
        }
        else {
            spec_values(k) = (std::complex<double>) 1j * (double) (k - n) * spec_values(k);
        }
    }
    fftw_execute(p_inverse_hanning);
    hanning_values.write("output/hann_" + std::to_string(n) + ".dat");

    return 0;
}

double f(double x) {
    return exp(sin(x));
}

double dfdx(double x) {
    return exp(sin(x)) * cos(x);
}

void fem_elem_matrix(double dx, double f1, double f2, Array<double>* A, Array<double>* b) {
    (*A)(0, 0) = dx / 3.0;
    (*A)(0, 1) = dx / 6.0;
    (*A)(1, 0) = dx / 6.0;
    (*A)(1, 1) = dx / 3.0;

    (*b)(0) = 1.0 / 2.0 * (f2 - f1);
    (*b)(1) = 1.0 / 2.0 * (f2 - f1);
}

void hanning_window(int N, Array<double>* w) {
    for (int i = 0; i < N; i++) {
        (*w)(i) = 0.5 * (1 + cos(2 * PI * (double) i / (double) N));
    }
}
