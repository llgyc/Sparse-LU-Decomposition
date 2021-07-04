#include <omp.h>
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
#include <unordered_map>

#define START(prompt) prof.starttime(prompt)
#define STOP(prompt) prof.stoptime(prompt)
#define THREAD_NUM 12
#define THREAD_NUM2 8
#define CHUNK_SIZE 36

int N, nonzero;
std::vector<std::vector<int>> row;
std::vector<std::vector<std::pair<int, double>>> col;

class DSU {
public:
    DSU(int size) { fa.resize(size); for (int i = 0; i < size; ++i) fa[i] = i; }
    int getfather(int x) { return (fa[x] == x) ? x : (fa[x] = getfather(fa[x])); }
    void merge(int x, int y) { fa[getfather(x)] = y; }
    bool same(int x, int y) { return getfather(x) == getfather(y); }

private:
    std::vector<int> fa;
};

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

std::vector<int> f;
// Column first
std::vector<std::queue<int>> tmp[THREAD_NUM2];
std::vector<std::vector<int>> id;
std::vector<std::vector<double>> mat;

void LU_Decomposition() {
    START("Elimination Tree");
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
    STOP("Elimination Tree");

    START("Fill-in");
    // Calculate fill-in
    id.resize(N); mat.resize(N);
    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) shared(tmp) num_threads(THREAD_NUM2)
    for (int i = 0; i < N; ++i) {
        std::vector<bool> tag; tag.resize(N);
        int myid = omp_get_thread_num();
        tmp[myid].resize(N);
        for (auto x : row[i]) {
            int now = x;
            while (now != i) {
                if (tag[now]) break;
                tag[now] = true;
                tmp[myid][now].push(i);
                now = f[now];
            }
        }
        tmp[myid][i].push(i);
    }
    row.clear(); f.clear();
    STOP("Fill-in");

    START("Initval");
    // Matrix Initialization
    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) shared(tmp, id, mat) num_threads(THREAD_NUM2)
    for (int i = 0; i < N; ++i) {
        unsigned ptr = 0;
        #define METHOD2
        #ifdef METHOD1
        std::priority_queue<std::pair<int, int>> aux;
        for (int j = 0; j < THREAD_NUM2; ++j) {
            if (tmp[j].size() == 0) continue;
            if (!tmp[j][i].empty()) {
                aux.push(std::make_pair(-tmp[j][i].front(), j));
                tmp[j][i].pop();
            }
        }
        while (!aux.empty()) {
            if (ptr < col[i].size() && col[i][ptr].first == -aux.top().first) {
                id[i].push_back(col[i][ptr].first);
                mat[i].push_back(col[i][ptr].second);
                ++ptr;
            } else {
                id[i].push_back(-aux.top().first);
                mat[i].push_back(0);
            }
            unsigned nowid = aux.top().second;
            if (!tmp[nowid][i].empty()) {
                aux.push(std::make_pair(-tmp[nowid][i].front(), nowid));
                tmp[nowid][i].pop();
            }
            aux.pop();
        }
        #endif
        #ifdef METHOD2
        while (true) {
            int pos = N, nowid = -1;
            for (int j = 0; j < THREAD_NUM2; ++j) {
                if (tmp[j].size() == 0) continue;
                if (!tmp[j][i].empty()) {
                    if (tmp[j][i].front() < pos) {
                        pos = tmp[j][i].front();
                        nowid = j;
                    }
                }
            }
            if (pos == N) break;
            if (ptr < col[i].size() && col[i][ptr].first == pos) {
                id[i].push_back(pos);
                mat[i].push_back(col[i][ptr].second);
                ++ptr;
            } else {
                id[i].push_back(pos);
                mat[i].push_back(0);
            }
            tmp[nowid][i].pop();
        }
        #endif
    }
    col.clear(); 
    for (int i = 0; i < THREAD_NUM2; ++i) tmp[i].clear();
    STOP("Initval");

    START("Bruteforce");
    // Bruteforce
    for (int k = 0; k < N; ++k) {
        double pivot = 1.0 / sqrt(mat[k][0]);
        int sz = mat[k].size();
        for (int j = 0; j < sz; ++j) {
            mat[k][j] *= pivot;
        }
        #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) shared(mat) num_threads(THREAD_NUM)
        for (int j = 1; j < sz; ++j) {
            double tmp = mat[k][j];
            int pos = id[k][j], idx = 0;
            for (int i = j; i < sz; ++i) {
                while (id[pos][idx] < id[k][i]) ++idx;
                mat[pos][idx] -= tmp * mat[k][i];
            }
        }
    }
    STOP("Bruteforce");

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
    row.resize(N); col.resize(N);
    for (int i = 0; i < nonzero; ++i) {
        int x, y;
        double val;
        fscanf(ifs, "%d%d%lf", &x, &y, &val);
        --x, --y;
        if (x != y) row[x].push_back(y);
        col[y].push_back(std::make_pair(x, val));
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
    for (int j = 0; j < N; ++j) num += id[j].size();
    //ofs_l << N << " " << N << " " << num << std::endl;
    fprintf(ofs_l, "%d %d %d\n", N, N, num);
    for (int j = 0; j < N; ++j) {
        if (id[j].size() == 0) continue;
        double pivot = mat[j][0];
        for (unsigned t = 0; t < id[j].size(); ++t) {
            int i = id[j][t];
            double val = mat[j][t] / pivot;
            //ofs_l << i + 1 << " " << j + 1 << " " << val << std::endl;
            fprintf(ofs_l, "%d %d %.15e\n", i + 1, j + 1, val);
        }
    }
    //ofs_l.close();
    fclose(ofs_l);

    //ofs_u << std::setprecision(16);
    //ofs_u << N << " " << N << " " << num << std::endl;
    fprintf(ofs_u, "%d %d %d\n", N, N, num);
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
            //ofs_u << i + 1 << " " << j + 1 << " " << val << std::endl;
            fprintf(ofs_u, "%d %d %.15e\n", i + 1, j + 1, val);
        }
    }
    //ofs_u.close();
    fclose(ofs_u);

    STOP("Total");
    
    return 0;
}
