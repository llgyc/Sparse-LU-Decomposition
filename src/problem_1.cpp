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
#include <unordered_map>

#define START(prompt) prof.starttime(prompt)
#define STOP(prompt) prof.stoptime(prompt)
class Profiling {
public:
    void starttime(const std::string &name) {
        start_time_[name] = std::chrono::steady_clock::now();
        fprintf(stderr, "[START] %s\n", name.c_str());
        //std::cerr << "[START] " << name << std::endl;
    }
    void stoptime(const std::string &name) {
        auto stop_time_ = std::chrono::steady_clock::now();
        /*std::cerr << "[STOP] " << name << " -- Time: " 
            << (std::chrono::duration_cast<std::chrono::microseconds>(stop_time_ - start_time_[name]).count()) 
            / 1000000.0 << "[s]" << std::endl;*/
        fprintf(stderr, "[STOP] %s -- Time: %.6lf[s]\n", name.c_str(), 
            (std::chrono::duration_cast<std::chrono::microseconds>(stop_time_ - start_time_[name]).count()) 
            / 1000000.0);
    }

private:
    std::unordered_map<std::string, std::chrono::steady_clock::time_point> start_time_;
} prof;

int N, nonzero;
/* Column-first */
std::vector<double> A;

/* Gauss Elimination */
void LU_Decomposition() {
    for (int k = 0; k < N - 1; ++k) {
        // Sanity check
        double pivot = A[k * N + k];
        if (fabs(pivot) == 0) {
            std::cerr << "[Erroe] Matrix A cannot be decomposed into form A = L * U" << std::endl;
            exit(-1);
        }
        for (int i = k + 1; i < N; ++i) {
            A[k * N + i] /= pivot;
        }
        for (int j = k + 1; j < N; ++j) {
            double tmp = A[j * N + k];
            for (int i = k + 1; i < N; ++i) {
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
    
    FILE *ifs = fopen(argv[1], "r");
    if (!ifs) {
        fprintf(stderr, "[ERROR] Cannot open matrix file\n");
        return -1;
    }
    /*
    std::ifstream ifs(matrix_path);
    if (!ifs.is_open()) {
        std::cerr << "[ERROR] Cannot open matrix file" << std::endl;
        return -1;
    }
    */

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

    START("Total");

    // Input
    /*std::string line;
    while (std::getline(ifs, line)) {
        if (line.size() > 0 && line[0] == '%')
            continue;
        std::istringstream iss(line);
        int M;
        iss >> N >> M >> nonzero;
        assert(N == M);
        break;
    }
    
    row.resize(N); col.resize(N);
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
        double val;
        ifs >> x >> y >> val;
        --x, --y;
        if (x != y) row[x].push_back(y);
        col[y].push_back(std::make_pair(x, val));
    }

    ifs.close();*/
    char ch = getc(ifs);
    while (ch == '%') {
        while (ch != '\n') ch = getc(ifs);
        ch = getc(ifs);
    }
    ungetc(ch, ifs);
    int M;
    fscanf(ifs, "%d%d%d", &N, &M, &nonzero);
    assert(N == M);
    A.resize(N * N);
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
        double val;
        fscanf(ifs, "%d%d%lf", &x, &y, &val);
        --x, --y;
        A[x * N + y] = val;
        A[y * N + x] = val;
    }
    fclose(ifs);

    // Calculation
    START("Calculation");
    LU_Decomposition();
    STOP("Calculation");
    
    // Output
    auto output_file1 = target_dir + "/2_L_" + matrix_name;
    auto output_file2 = target_dir + "/2_U_" + matrix_name;
    FILE *ofs_l = fopen(output_file1.c_str(), "w");
    FILE *ofs_u = fopen(output_file2.c_str(), "w");
    //std::ofstream ofs_l(output_file1);
    //std::ofstream ofs_u(output_file2);

    //ofs_l << std::setprecision(16);
    int num = 0;
    for (int j = 0; j < N; ++j) {
        for (int i = j; i < N; ++i) {
            if (i == j || fabs(A[j * N + i]) != 0) {
                ++num;
            }
        }
    }
    fprintf(ofs_l, "%d %d %d\n", N, N, num);
    for (int j = 0; j < N; ++j) {
        for (int i = j; i < N; ++i) {
            if (i == j) 
                fprintf(ofs_l, "%d %d 1\n", i + 1, j + 1);
            else if (fabs(A[j * N + i]) != 0)
                fprintf(ofs_l, "%d %d %.15e\n", i + 1, j + 1, A[j * N + i]);
        }
    }
    //ofs_l.close();
    fclose(ofs_l);

    //ofs_u << std::setprecision(16);
    //ofs_u << N << " " << N << " " << num << std::endl;
    num = 0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i <= j; ++i) {
            if (fabs(A[j * N + i]) != 0) {
                ++num;
            }
        }
    }
    fprintf(ofs_u, "%d %d %d\n", N, N, num);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i <= j; ++i) {
            if (fabs(A[j * N + i]) != 0)
                fprintf(ofs_u, "%d %d %.15e\n", i + 1, j + 1, A[j * N + i]);
        }
    }    
    //ofs_u.close();
    fclose(ofs_u);
    
    STOP("Total");

    return 0;
}
