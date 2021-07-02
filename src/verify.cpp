#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

const double eps = 1e-6;

int size = -1;
std::vector<double> A, L, U, result;

void read_matrix(std::ifstream &ifs, std::vector<double> &mat, bool copy = false) {
    int N, nonzero;
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
    if (copy) size = N;
    else assert(N == size);
    mat.resize(N * N);
    for (int i = 0; i < N * N; ++i) mat[i] = 0;
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
        double val;
        ifs >> x >> y >> val;
        --x, --y;
        mat[x * N + y] = val;
        if (copy && x != y) mat[y * N + x] = val;
    }
    ifs.close();
}

void multiply(std::vector<double> &C, std::vector<double> &A, std::vector<double> &B) {
    std::vector<std::vector<int>> tmp;
    tmp.resize(size);
    for (int k = 0; k < size; ++k) {
        tmp[k] = std::vector<int>();
        for (int j = 0; j < size; ++j) {
            if (B[k * size + j] != 0) tmp[k].push_back(j);
        }
    }

    for (int i = 0; i < size * size; ++i) C[i] = 0;
    for (int i = 0; i < size; ++i) {
        for (int k = 0; k < size; ++k) {
            if (A[i * size + k] == 0) continue;
            for (auto x : tmp[k]) {
                C[i * size + x] += A[i * size + k] * B[k * size + x];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    // Environment setup
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [MATRIX A] [MATRIX L] [MATRIX U]" << std::endl;
        return -1;
    }
    std::ifstream ifs_a(argv[1]);
    read_matrix(ifs_a, A, true);

    std::ifstream ifs_l(argv[2]);
    read_matrix(ifs_l, L);

    std::ifstream ifs_u(argv[3]);
    read_matrix(ifs_u, U);
    
    result.resize(size * size);
    multiply(result, L, U);

    double mx = 0;
    for (int q = 0; q < size; ++q) {
        for (int p = q; p < size; ++p) {
            int i = p * size + q;
            double tmp = fabs(result[i]);
            tmp = std::min(tmp, fabs(result[i] - A[i]) / A[i]);
            if (tmp > 1e-5)
                std::cerr << i << ": " << result[i] << " " << A[i] << std::endl;
            mx = std::max(mx, tmp);
        }
    }
    std::cout << "Max Relative Error (Single point): " << mx << std::endl;

    return 0;
}
