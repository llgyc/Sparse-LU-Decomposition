#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <queue>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

int N;
int nonzero;
std::vector<std::vector<int>> row;
std::vector<std::vector<std::pair<int, double>>> col;

class DSU {
public:
    DSU(int size) { fa.resize(size); for (int i = 0; i < size; ++i) fa[i] = i; }
    int getfather(int x) { return (fa[x] == x) ? x : (fa[x] = getfather(fa[x])); }
    void merge(int x, int y) { fa[getfather(x)] = y; }
    bool same(long x, long y) { return getfather(x) == getfather(y); }

private:
    std::vector<int> fa;
};

std::vector<int> f, tag;
// Column first
std::vector<std::vector<int>> id;
std::vector<std::vector<double>> mat;

void LU_Decomposition() {
    // Build Elimination Tree
    DSU dsu(N); f.resize(N);
    for (auto &x : f) x = -1;
    for (int i = 0; i < N; ++i) {
        for (auto j : row[i]) {
            if (dsu.same(i,j)) continue;
            int root = dsu.getfather(j);
            f[root] = i;
            dsu.merge(root, i);
        }
    }

    // Calculate fill-in
    tag.resize(N); id.resize(N); mat.resize(N);
    for (auto &x : tag) x = -1;
    for (int i = 0; i < N; ++i) {
        for (auto x : row[i]) {
            int now = x;
            while (now != i) {
                if (tag[now] == i) break;
                tag[now] = i;
                id[now].push_back(i);
                mat[now].push_back(0);
                now = f[now];
            }
        }
        id[i].push_back(i);
        mat[i].push_back(0);
    }
    row.clear(); tag.clear(); f.clear(); 

    // Matrix Initialization
    for (int i = 0; i < N; ++i) {
        unsigned ptr = 0;
        for (unsigned j = 0; j < id[i].size(); ++j) {
            if (ptr < col[i].size()) { 
                if (col[i][ptr].first == id[i][j]) {
                    mat[i][j] = col[i][ptr].second; ++ptr;
                
                }
            } else break;
        }
    }
    col.clear(); 

    // Bruteforce
    for (int k = 0; k < N; ++k) {
        double pivot = 1.0 / sqrt(mat[k][0]);
        int sz = mat[k].size();
        for (int j = 0; j < sz; ++j) {
            mat[k][j] *= pivot;
        }
        for (int j = 1; j < sz; ++j) {
            double tmp = mat[k][j];
            int pos = id[k][j], idx = 0;
            for (int i = j; i < sz; ++i) {
                while (id[pos][idx] < id[k][i]) ++idx;
                mat[pos][idx] -= tmp * mat[k][i];
            }
        }
    }

}

std::vector<std::vector<std::pair<int, double>>> out;

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
    
    row.resize(N); col.resize(N);
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
        double val;
        ifs >> x >> y >> val;
        --x, --y;
        if (x != y) row[x].push_back(y);
        col[y].push_back(std::make_pair(x, val));
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
    /*
    auto output_file1 = target_dir + "/2_L_" + matrix_name;
    auto output_file2 = target_dir + "/2_U_" + matrix_name;
    std::ofstream ofs_l(output_file1);
    std::ofstream ofs_u(output_file2);

    ofs_l << std::setprecision(16);
    long long num = 0;
    for (int j = 0; j < N; ++j) num += id[j].size();
    ofs_l << N << " " << N << " " << num << std::endl;
    for (int j = 0; j < N; ++j) {
        if (id[j].size() == 0) continue;
        double pivot = mat[j][0];
        for (unsigned t = 0; t < id[j].size(); ++t) {
            int i = id[j][t];
            double val = mat[j][t] / pivot;
            ofs_l << i + 1 << " " << j + 1 << " " << val << std::endl;
        }
    }
    ofs_l.close();

    ofs_u << std::setprecision(16);
    ofs_u << N << " " << N << " " << num << std::endl;
    out.resize(N);
    for (int j = 0; j < N; ++j) {
        if (id[j].size() == 0) continue;
        double pivot = mat[j][0];
        for (unsigned t = 0; t < id[j].size(); ++t) {
            out[id[j][t]].push_back(std::make_pair(j, mat[j][t] * pivot));
        }
    }
    for (int j = 0; j < N; ++j) {
        for (auto x : out[j]) {
            int i = x.first;
            double val = x.second;
            ofs_u << i + 1 << " " << j + 1 << " " << val << std::endl;
        }
    }
    ofs_u.close();
    */
    auto finish = std::chrono::steady_clock::now();
    std::cerr << "Total Time: " 
        << (std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()) 
        / 1000000.0 << "[s]" << std::endl;
    
    return 0;
}
