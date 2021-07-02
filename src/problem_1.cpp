#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

long long N, nonzero;
/* Column-first */
std::vector<double> A;

/* Gauss Elimination */
void LU_Decomposition() {
    for (long long k = 0; k < N - 1; ++k) {
        // Sanity check
        double pivot = A[k * N + k];
        if (fabs(pivot) == 0) {
            std::cerr << "[Erroe] Matrix A cannot be decomposed into form A = L * U" << std::endl;
            exit(-1);
        }
        for (long long i = k + 1; i < N; ++i) {
            A[k * N + i] /= pivot;
        }
        for (long long j = k + 1; j < N; ++j) {
            double tmp = A[j * N + k];
            for (long long i = k + 1; i < N; ++i) {
                A[j * N + i] -= A[k * N + i] * tmp;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    // Environment setup
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [MATRIX FILE]" << std::endl;
        return -1;
    }
    
    std::string file_path(argv[0]);
    std::string matrix_path(argv[1]);
    
    std::ifstream ifs(matrix_path);
    if (!ifs.is_open()) {
        std::cerr << "[ERROR] Cannot open matrix file" << std::endl;
        return -1;
    }

    auto file_pos = file_path.find_last_of("/");
    std::string target_dir;
    if (file_pos == std::string::npos) {
        target_dir = "../result";
    } else {
        target_dir = file_path.substr(0, file_pos) + "/../result";
    }

    auto matrix_pos = matrix_path.find_last_of("/");
    std::string matrix_name = (matrix_pos == std::string::npos) 
        ? matrix_path : matrix_path.substr(matrix_pos + 1);

    struct stat st = {0};
    if (stat(target_dir.c_str(), &st) == -1) {
        mkdir(target_dir.c_str(), 0700);
    }

    auto start = std::chrono::steady_clock::now();

    // Input
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.size() > 0 && line[0] == '%')
            continue;
        std::istringstream iss(line);
        long long M;
        iss >> N >> M >> nonzero;
        assert(N == M);
        break;
    }
    A.resize(N * N);
    for (long long i = 0; i < N * N; ++i) A[i] = 0;
    for (long long i = 0; i < nonzero; ++i) {
        long long x, y;
        double val;
        ifs >> x >> y >> val;
        --x, --y;
        A[x * N + y] = val;
        A[y * N + x] = val;
    }
    ifs.close();

    // Calculation
    auto begin  = std::chrono::steady_clock::now();

    LU_Decomposition();

    auto end = std::chrono::steady_clock::now();
    std::cerr << "Calculation Time: " 
        << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) 
        / 1000000.0 << "[s]" << std::endl;

    // Output
    auto output_file1 = target_dir + "/1_L_" + matrix_name;
    auto output_file2 = target_dir + "/1_U_" + matrix_name;
    std::ofstream ofs_l(output_file1);
    std::ofstream ofs_u(output_file2);

    ofs_l << std::setprecision(16);
    long long num = 0;
    for (long long j = 0; j < N; ++j) {
        for (long long i = j; i < N; ++i) {
            if (i == j || fabs(A[j * N + i]) != 0) {
                ++num;
            }
        }
    }
    ofs_l << N << " " << N << " " << num << std::endl;
    for (long long j = 0; j < N; ++j) {
        for (long long i = j; i < N; ++i) {
            if (i == j)
                ofs_l << i + 1 << " " << j + 1 << " 1" << std::endl;
            else if (fabs(A[j * N + i]) != 0)
                ofs_l << i + 1 << " " << j + 1 << " " << A[j * N + i] << std::endl;
        }
    }
    ofs_l.close();

    ofs_u << std::setprecision(16);
    num = 0;
    for (long long j = 0; j < N; ++j) {
        for (long long i = 0; i <= j; ++i) {
            if (fabs(A[j * N + i]) != 0) {
                ++num;
            }
        }
    }
    ofs_u << N << " " << N << " " << num << std::endl;
    for (long long j = 0; j < N; ++j) {
        for (long long i = 0; i <= j; ++i) {
            if (fabs(A[j * N + i]) != 0)
                ofs_u << i + 1 << " " << j + 1 << " " << A[j * N + i] << std::endl;
        }
    }
    ofs_u.close();

    auto finish = std::chrono::steady_clock::now();
    std::cerr << "Total Time: " 
        << (std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()) 
        / 1000000.0 << "[s]" << std::endl;

    return 0;
}
