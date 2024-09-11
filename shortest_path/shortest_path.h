#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <set>

class Graph {
public:
    long long cnt=0;
    int number_of_nodes;  
    std::vector<double> lon, lat;
    int number_of_line;
    std::vector<std::vector<std::pair<int, int>>> graph;

    std::vector<int>pre;
    std::vector<int>distance;
    std::vector<int>predict;
    std::vector<int>cnt_of_path;
    std::map<int, std::map<int, bool>>exist_line;

    Graph(const std::string& node_file, const std::string& edge_file);
    
    void readin(const std::string& node_file, const std::string& edge_file);
   
    void self_test();
    void check();
    int number_of_query;
    void work();
    int ALT(int s,int t);
    int Astar(int s,int t);
    int dijkstra(int s,int t);
    std::vector<int>frontdis,backdis;
    int Bidijkstra(int s,int t);
    double degreesToRadians(double degrees);
    int real_dis(double lat1, double lon1, double lat2, double lon2);
    int query_h(int node,int t);
    int CH(int s,int t);

    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> dijkstra_way(int s,int t);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> bi_dijkstra_way(int s,int t);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> ALT_way(int s,int t);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Astar_way(int s,int t);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> CH_way(int s,int t);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> H2H_way(int s,int t,int dist);
    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> PLL_way(int s,int t,int dist);
    struct result {
        int length,s,t;
    };
    static bool compare(const result& a, const result& b) {
        return a.length < b.length;
    }
    std::vector<result>test;

    std::vector<std::vector<int>> Preprocessed_path;
    void preprocess(int order,int s);
    std::vector<std::vector<std::pair<int,int>>>queries;
    void makequery();

    struct LabelEntry {
        int vertex;
        int distance;
    };
    struct node{
        int osmid,degree;
    };
    static bool comparenodes(const node& a, const node& b);
    static bool compareLabel(const LabelEntry& a, const LabelEntry& b);
    std::vector<node> nodes;
    std::vector<int> order;
    std::vector<std::vector<LabelEntry>> Labels;
    void preprocess2Hop();
    int q2Hop(int u, int v);
    void prunedBFS(int start);

    //CH
    void CH_preprocess();
    int mx_dis,sec_mx;
    std::vector<int>CH_order;
    std::vector<int>hf_dis;
    std::vector<int> contracted_nodes;
    std::map<std::pair<int, int>, std::vector<int> > shortcut_map;

    int tim;
    std::vector<int>time_ch;
    int compute_edge_difference(int node);
    int compute_deleted_neighbors(int node);

    //H2H
    void H2H_preprocess();
    int lgsiz;
    std::vector<std::vector<int>>dis;
    std::vector<std::vector<int>>pos;
    std::vector<std::vector<int>>X;
    std::vector<std::vector<int>>fi_dis;
    std::vector<std::vector<int>>ancestor;
    std::vector<int>H2H_order;
    std::vector<std::vector<int>>fa;
    std::vector<int>dep;
    std::vector<std::vector<std::pair<int, int>>> H2H_graph;
    int H2H(int s,int t);
    int LCA(int s,int t);
    bool compareX(int a, int b);
    std::vector<int> H2H_path(int s,int t,int dis);

    void save_preprocessed_data(
        const std::vector<int>& CH_order, 
        const std::map<std::pair<int, int>, std::vector<int>>& shortcut_map,
        const std::vector<int>& hf_dis,
        const std::vector<int>& contracted_nodes) {

        std::ofstream outfile("preprocessed_data.bin", std::ios::binary);

        // Save CH_order
        size_t order_size = CH_order.size();
        outfile.write(reinterpret_cast<const char*>(&order_size), sizeof(order_size));
        outfile.write(reinterpret_cast<const char*>(CH_order.data()), order_size * sizeof(int));

        // Save shortcut_map
        size_t map_size = shortcut_map.size();
        outfile.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        for (const auto& entry : shortcut_map) {
            outfile.write(reinterpret_cast<const char*>(&entry.first.first), sizeof(int));
            outfile.write(reinterpret_cast<const char*>(&entry.first.second), sizeof(int));
            
            size_t vector_size = entry.second.size();
            outfile.write(reinterpret_cast<const char*>(&vector_size), sizeof(vector_size));
            outfile.write(reinterpret_cast<const char*>(entry.second.data()), vector_size * sizeof(int));
        }

        // Save hf_dis
        size_t hf_dis_size = hf_dis.size();
        outfile.write(reinterpret_cast<const char*>(&hf_dis_size), sizeof(hf_dis_size));
        outfile.write(reinterpret_cast<const char*>(hf_dis.data()), hf_dis_size * sizeof(int));

        // Save contracted_nodes
        size_t contracted_nodes_size = contracted_nodes.size();
        outfile.write(reinterpret_cast<const char*>(&contracted_nodes_size), sizeof(contracted_nodes_size));
        outfile.write(reinterpret_cast<const char*>(contracted_nodes.data()), contracted_nodes_size * sizeof(int));

        outfile.close();
    }

    bool load_preprocessed_data(
        std::vector<int>& CH_order, 
        std::map<std::pair<int, int>, std::vector<int>>& shortcut_map,
        std::vector<int>& hf_dis,
        std::vector<int>& contracted_nodes) {

        std::ifstream infile("preprocessed_data.bin", std::ios::binary);
        if (!infile) return false;
        
        // Load CH_order
        size_t order_size;
        infile.read(reinterpret_cast<char*>(&order_size), sizeof(order_size));
        CH_order.resize(order_size);
        infile.read(reinterpret_cast<char*>(CH_order.data()), order_size * sizeof(int));

        // Load shortcut_map
        size_t map_size;
        infile.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        for (size_t i = 0; i < map_size; ++i) {
            int first, second;
            infile.read(reinterpret_cast<char*>(&first), sizeof(int));
            infile.read(reinterpret_cast<char*>(&second), sizeof(int));

            size_t vector_size;
            infile.read(reinterpret_cast<char*>(&vector_size), sizeof(vector_size));
            std::vector<int> shortcut(vector_size);
            infile.read(reinterpret_cast<char*>(shortcut.data()), vector_size * sizeof(int));

            shortcut_map[{first, second}] = shortcut;
        }

        // Load hf_dis
        size_t hf_dis_size;
        infile.read(reinterpret_cast<char*>(&hf_dis_size), sizeof(hf_dis_size));
        hf_dis.resize(hf_dis_size);
        infile.read(reinterpret_cast<char*>(hf_dis.data()), hf_dis_size * sizeof(int));

        // Load contracted_nodes
        size_t contracted_nodes_size;
        infile.read(reinterpret_cast<char*>(&contracted_nodes_size), sizeof(contracted_nodes_size));
        contracted_nodes.resize(contracted_nodes_size);
        infile.read(reinterpret_cast<char*>(contracted_nodes.data()), contracted_nodes_size * sizeof(int));

        infile.close();
        return true;
    }


};

#define NULLINDEX 0xFFFFFFFF
class heap {
public:
    struct element_t {
        int key;
        int element;
        element_t() : key(0), element(0) {}
        element_t(const int k, const int e) : key(k), element(e) {}
    };

private:
    int n;
    int max_n;
    std::vector<element_t> elements;
    std::vector<int> position;

public:
    heap(int n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {}
    heap() {}
    inline int size() const { return n; }
    inline bool empty() const { return size() == 0; }
    inline void extract_min(int &element, int &key) {
        assert(!empty());
        element_t &front = elements[0];
        element = front.element;
        key = front.key;
        position[element] = NULLINDEX;
        --n;
        if (!empty()) {
            front = elements[n];
            position[front.element] = 0;
            sift_down(0);
        }
    }
    inline int top() {
        assert(!empty());
        return elements[0].key;
    }
    inline int top_value() {
        assert(!empty());
        return elements[0].element;
    }
    inline void update(const int element, const int key) {
        if (position[element] == NULLINDEX) {
            element_t &back = elements[n];
            back.key = key;
            back.element = element;
            position[element] = n;
            sift_up(n++);
        } else {
            int el_pos = position[element];
            element_t &el = elements[el_pos];
            if (key > el.key) {
                el.key = key;
                sift_down(el_pos);
            } else {
                el.key = key;
                sift_up(el_pos);
            }
        }
    }
    inline void clear() {
        for (int i = 0; i < n; ++i) {
            position[elements[i].element] = NULLINDEX;
        }
        n = 0;
    }
    inline bool contains(const int element) const {
        return position[element] != NULLINDEX;
    }

protected:
    inline void sift_up(int i) {
        assert(i < n);
        int cur_i = i;
        while (cur_i > 0) {
            int parent_i = (cur_i-1) >> 1;
            if (elements[parent_i].key > elements[cur_i].key) {
                swap(cur_i, parent_i);
            } else {
                break;
            }
            cur_i = parent_i;
        }
    }
    inline void sift_down(int i) {
        assert(i < n);
        while (true) {
            int min_ind = i;
            int min_key = elements[i].key;
            int child_ind_l = (i << 1) + 1;
            int child_ind_u = std::min(child_ind_l + 2, n);
            for (int j = child_ind_l; j < child_ind_u; ++j) {
                if (elements[j].key < min_key) {
                    min_ind = j;
                    min_key = elements[j].key;
                }
            }
            if (min_ind != i) {
                swap(i, min_ind);
                i = min_ind;
            } else {
                break;
            }
        }
    }
    inline void swap(const int i, const int j) {
        element_t &el_i = elements[i];
        element_t &el_j = elements[j];
        position[el_i.element] = j;
        position[el_j.element] = i;
        element_t temp = el_i;
        el_i = el_j;
        el_j = temp;
    }
};
#endif //SHORTEST_PATH_H
