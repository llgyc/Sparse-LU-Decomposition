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

long long N, nonzero;
std::vector<std::vector<long long>> row;
std::vector<std::vector<std::pair<long long, double>>> col;

class DSU {
public:
    DSU(long long size) { fa.resize(size); for (long long i = 0; i < size; ++i) fa[i] = i; }
    long long getfather(long long x) { return (fa[x] == x) ? x : (fa[x] = getfather(fa[x])); }
    void merge(long long x, long long y) { fa[getfather(x)] = y; }
    bool same(long x, long y) { return getfather(x) == getfather(y); }

private:
    std::vector<long long> fa;
};

std::vector<long long> f, tag;
std::vector<std::vector<long long>> edges;
// Column first
std::vector<std::vector<double>> id;
std::vector<std::vector<double>> mat;

void CalculateDense(long long x) {
    long long sz = id[x].size();
    double pivot = 1.0 / sqrt(mat[x][0]);
    for (long long i = 0; i < sz; ++i) {
        mat[x][i] *= pivot;
    }
    for (long long j = 1; j < sz; ++j) {
        double tmp = mat[x][j];
        for (long long i = j; i < sz; ++i) {
            mat[x][j * sz + i] -= mat[x][i] * tmp;
        }
    }
}

void UpdateMatrix(long long x) {
    std::vector<long long> rel;
    long long fa = f[x];
    unsigned ptr = 1;
    for (unsigned i = 0; i < id[fa].size(); ++i) {
        if (ptr < id[x].size()) {
            if (id[fa][i] == id[x][ptr]) {
                rel.push_back(i);
                ++ptr;
            }
        } else break;
    }
    
    long long sz = id[x].size();
    long long sz2 = id[fa].size();
    for (long long j = 1; j < sz; ++j) {
        for (long long i = j; i < sz; ++i) {
            mat[fa][rel[j-1] * sz2 + rel[i-1]] += mat[x][j * sz + i];
        }
    }
}

void Multifrontal(long long x) {
    mat[x].resize(id[x].size() * id[x].size());
    // Postorder
    for (auto y : edges[x]) {
        Multifrontal(y);
    }
    CalculateDense(x);
    UpdateMatrix(x);
    mat[x].resize(id[x].size());
}

void LU_Decomposition() {
    // Build Elimination Tree
    DSU dsu(N); f.resize(N); edges.resize(N);
    for (auto &x : f) x = -1;
    for (long long i = 0; i < N; ++i) {
        for (auto j : row[i]) {
            if (dsu.same(i,j)) continue;
            long long root = dsu.getfather(j);
            f[root] = i;
            edges[i].push_back(root);
            dsu.merge(root, i);
        }
    }

    // Calculate fill-in
    tag.resize(N); id.resize(N); mat.resize(N);
    for (auto &x : tag) x = -1;
    for (long long i = 0; i < N; ++i) {
        for (auto x : row[i]) {
            long long now = x;
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
    row.clear();

    // Matrix Initialization
    for (long long i = 0; i < N; ++i) {
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
    
    // Multifrontal
    for (long long i = 0; i < N; ++i)
        if (f[i] == -1) Multifrontal(i);

}

std::vector<std::vector<std::pair<long long, double>>> out;

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
    
    row.resize(N); col.resize(N);
    for (long long i = 0; i < nonzero; ++i) {
        long long x, y;
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
    auto output_file1 = target_dir + "/2_L_" + matrix_name;
    auto output_file2 = target_dir + "/2_U_" + matrix_name;
    std::ofstream ofs_l(output_file1);
    std::ofstream ofs_u(output_file2);

    // TODO
    ofs_l << std::setprecision(16);
    long long num = 0;
    for (long long j = 0; j < N; ++j) num += id[j].size();
    ofs_l << N << " " << N << " " << num << std::endl;
    for (long long j = 0; j < N; ++j) {
        if (id[j].size() == 0) continue;
        double pivot = mat[j][0];
        for (unsigned t = 0; t < id[j].size(); ++t) {
            long long i = id[j][t];
            double val = mat[j][t] / pivot;
            ofs_l << i + 1 << " " << j + 1 << " " << val << std::endl;
        }
    }
    ofs_l.close();

    ofs_u << std::setprecision(16);
    ofs_u << N << " " << N << " " << num << std::endl;
    out.resize(N);
    for (long long j = 0; j < N; ++j) {
        if (id[j].size() == 0) continue;
        double pivot = mat[j][0];
        for (unsigned t = 0; t < id[j].size(); ++t) {
            out[id[j][t]].push_back(std::make_pair(j, mat[j][t] * pivot));
        }
    }
    for (long long j = 0; j < N; ++j) {
        for (auto x : out[j]) {
            long long i = x.first;
            double val = x.second;
            ofs_u << i + 1 << " " << j + 1 << " " << val << std::endl;
        }
    }
    ofs_u.close();

    auto finish = std::chrono::steady_clock::now();
    std::cerr << "Total Time: " 
        << (std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()) 
        / 1000000.0 << "[s]" << std::endl;

    return 0;
}
