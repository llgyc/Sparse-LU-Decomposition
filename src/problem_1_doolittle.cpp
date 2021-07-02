#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

const double eps = 1e-6;

int N, nonzero;
std::vector<double> A, L, U;

void LU_Decomposition() {
    for (int i = 0; i < N; ++i) {
        // Calculate U first
        for (int j = i; j < N; ++j) {
            double temp = 0;
            for (int k = 0; k < i - 1; ++k) {
                temp += L[i * N + k] * U[k * N + j];
            }
            U[i * N + j] = A[i * N + j] - temp;
        }

        // Sanity check
        if (fabs(U[i * N + i]) <= eps) {
            std::cerr << "[Erroe] Matrix A cannot be decomposed into form A = L * U" << std::endl;
            exit(-1);
        }

        // Calculate L next
        for (int j = i + 1; j < N; ++j) {
            double temp = 0;
            for (int k = 0; k < i - 1; ++k) {
                temp += L[j * N + k] * U[k * N + i];
            }
            L[j * N + i] = (A[j * N + i] - temp) / U[i * N + i];
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
        int M;
        iss >> N >> M >> nonzero;
        assert(N == M);
        break;
    }
    A.resize(N * N);
    L.resize(N * N);
    U.resize(N * N);
    for (int i = 0; i < N * N; ++i) A[i] = L[i] = U[i] = 0;
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
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

    int num = 0;
    for (int i = 0; i < N * N; ++i)
        if (fabs(L[i]) > eps) ++num;
    ofs_l << N << " " << N << " " << num << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (fabs(L[i * N + j]) > eps)
                ofs_l << i + 1 << " " << j + 1 << " " << L[i * N + j] << std::endl;
        }
    }

    ofs_l.close();

    num = 0;
    for (int i = 0; i < N * N; ++i)
        if (fabs(U[i]) > eps) ++num;
    ofs_u << N << " " << N << " " << num << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (fabs(U[i * N + j]) > eps)
                ofs_u << i + 1 << " " << j + 1 << " " << U[i * N + j] << std::endl;
        }
    }

    ofs_u.close();

    auto finish = std::chrono::steady_clock::now();
    std::cerr << "Total Time: " 
        << (std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()) 
        / 1000000.0 << "[s]" << std::endl;

    return 0;
}
