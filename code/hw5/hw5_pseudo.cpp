#include "hw5.hpp"
#define K 4.0
#define c 4.0
#define v 0.2
#define PI 3.1415926535897932384626433
#define TF PI / 4.0

int main(int argc, char** argv) {
    int N = 8;
    if (argc == 2) {
        N = std::stoi(argv[1]);
    }

    double dt = 8.0 / v / pow(N, 2);
    double m = TF / dt;
    int timesteps = ceil(m) * 20;
    dt = TF / (double) timesteps;

    std::vector<std::complex<double>> u(N, 0.0);
    std::vector<std::complex<double>> u_prime(N, 0.0);
    std::vector<std::complex<double>> u_dbl_prime(N, 0.0);

    std::vector<std::complex<double>> u_tilde(N, 0.0);
    std::vector<std::complex<double>> u_prime_tilde(N, 0.0);
    std::vector<std::complex<double>> u_dbl_prime_tilde(N, 0.0);

    std::vector<double> u_out(N, 0.0);

    // Create FFTW plans
    fftw_plan u_forward = fftw_plan_dft_1d(N, (fftw_complex*) u.data(), (fftw_complex*) u_tilde.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan u_prime_forward = fftw_plan_dft_1d(N, (fftw_complex*) u_prime.data(), (fftw_complex*) u_prime_tilde.data(), FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan u_inverse = fftw_plan_dft_1d(N, (fftw_complex*) u_tilde.data(), (fftw_complex*) u.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan u_prime_inverse = fftw_plan_dft_1d(N, (fftw_complex*) u_prime_tilde.data(), (fftw_complex*) u_prime.data(), FFTW_BACKWARD, FFTW_ESTIMATE);

    // Read ICs
    std::string filename = "problem_definition/u0_" + std::to_string(N) + ".dat";
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::invalid_argument("File (" + filename + ") could not be opened, check to make sure its spelled correctly");
    }

    for (int i = 0; i < N; i++) {
        std::string line;
        std::string token;
        getline(file, line);

        u[i] = std::stod(line);
    }

    file.close();

    // Pseudospectral
    fftw_execute(u_forward);

    for (int l = 0; l < timesteps; l++) {
        // Calculate u' and u'' in frequency domain
        for (int k = 0; k < N; k++) {
            if (k < N / 2) {
                u_prime_tilde[k] = std::complex<double>(0.0, 1.0) * (double) k * u_tilde[k];
                u_dbl_prime_tilde[k] = -(double) pow(k, 2) * u_tilde[k];
            }
            else {
                std::complex<double> coeff(0.0, (double) (k - N));
                u_prime_tilde[k] = coeff * u_tilde[k];
                u_dbl_prime_tilde[k] = pow(coeff, 2) * u_tilde[k];
            }
        }

        // Get u * u' in real domain
        fftw_execute(u_inverse);
        fftw_execute(u_prime_inverse);

        for (int j = 0; j < N; j++) {
            u[j] = u[j] / (double) N;
            u_prime[j] = u_prime[j] / (double) N;
        }

        // Calculate u * u'
        for (int i = 0; i < N; i++) {
            u_prime[i] = u_prime[i] * u[i];
        }

        // Get u * u' in frequency domain
        fftw_execute(u_prime_forward);

        // Step time
        for (int k = 0; k < N; k++) {
            u_tilde[k] = u_tilde[k] + dt * (v * u_dbl_prime_tilde[k] - u_prime_tilde[k]);
        }
    }

    fftw_execute(u_inverse);

    for (int j = 0; j < N; j++) {
        u_out[j] = u[j].real() / (double) N;
    }

    write_real_vector(N, u_out.data(), "output/uf_pseudo_pi_4_" + std::to_string(N) + ".dat");

}

void print_real_vector(int n, double* vec) {
    std::cout << "<";
    for (int i = 0; i < n-1; i++) {
        std::cout << std::to_string(vec[i]) << ", ";
    }
    std::cout << std::to_string(vec[n-1]) << ">" << std::endl;
}

void write_real_vector(int n, double* vec, std::string filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::invalid_argument("File (" + filename + ") could not be opened for output, check to make sure its spelled correctly");
    }
    
    for (int i = 0; i < n-1; i++) {
        file << vec[i] << std::endl;
    }
    file << vec[n-1];

    file.close();
}

void print_complex_vector(int n, std::complex<double>* vec) {
    std::cout << "<";
    for (int i = 0; i < n-1; i++) {
        std::cout << std::to_string(vec[i].real()) << " + i" << std::to_string(vec[i].imag()) << ", ";
    }
    std::cout << std::to_string(vec[n-1].real()) << " + i" << std::to_string(vec[n-1].imag()) << ">" << std::endl;
}
