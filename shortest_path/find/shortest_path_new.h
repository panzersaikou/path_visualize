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

    Graph();
    
    void readin();
    void self_test();
    void check();
    int number_of_query;
    void work();
    int ALT(int s,int t);
    int Astar(int s,int t);
    int dijkstra(int s,int t);
    std::vector<int> dijkstra_path(int s,int t);
    std::vector<int>frontdis,backdis;
    int Bidijkstra(int s,int t);
    double degreesToRadians(double degrees);
    int real_dis(double lat1, double lon1, double lat2, double lon2);
    int query_h(int node,int t);
    int CH(int s,int t);

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

    //PLL
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
    std::vector<int> PLL_path(int s,int t,int dis);

    //CH
    void CH_preprocess();
    int mx_dis,sec_mx;
    std::vector<int>CH_order;
    std::vector<int>hf_dis;
    std::vector<int> contracted_nodes;
    std::vector<std::vector<std::pair<int, int>>> CH_graph;
    int tim;
    std::vector<int>time_ch;

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

    //CRP
    void CRP_preprocess();
    int CRP(int s,int t);
    int CRP_dijkstra(int s,int t);
    void random_test();
    std::vector<std::vector<int> >partition;
    std::vector<int>belong;
    std::vector<int>is_border;
    std::vector<std::vector<std::pair<int, int>>> CRP_graph;
    std::vector<std::vector<std::pair<int, int>>> CRP_distance;
    int flag;

    //CCH
    void CCH_preprocess();
    int top_of_union(int u);
    void customization();
    int CCH(int s,int t);
    std::vector<int>union_find;
    std::vector<int>is_super;
    std::vector<int>deleted;
    std::vector<int>CCH_order;
    std::vector<int>CCH_nodes;
    std::vector<std::unordered_map<int,int>>CCH_times;

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
