#include "shortest_path.h"
#include <cstdio>
#include <queue>
#include <utility>
#include <limits>
#include <vector>
#include <stack>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <set>

const double PI = 3.14159265358979323846;

// Convert an Angle to a function of radians
double Graph::degreesToRadians(double degrees) {
    return degrees * PI / 180.0;
}

Graph::Graph(const std::string& node_file,const std::string& graph_file) 
{
    readin(node_file,graph_file);
}


bool Graph::comparenodes(const node& a, const node& b)
{
    return a.degree>b.degree;
}

bool Graph::compareLabel(const LabelEntry& a, const LabelEntry& b)
{
    return a.vertex<b.vertex;
}

bool Graph::compareX(int a, int b)
{
    return H2H_order[a]>H2H_order[b];
}

void Graph::readin(const std::string& node_file,const std::string& graph_file) 
{
    int bo=1;
    exist_line.clear();
    //readin a traffic graph
    std::ifstream file(node_file);
    long long ix,iy;
    double dx, dy;
    //get the information of nodes
    if (file.is_open()) {
        file >> number_of_nodes;
        lon.resize(number_of_nodes);
        lat.resize(number_of_nodes);
        pre.resize(number_of_nodes);
        distance.resize(number_of_nodes);
        graph.resize(number_of_nodes);
        predict.resize(number_of_nodes);
        frontdis.resize(number_of_nodes);
        backdis.resize(number_of_nodes);
        nodes.resize(number_of_nodes);
        Labels.resize(number_of_nodes);
        order.resize(number_of_nodes);
        contracted_nodes.resize(number_of_nodes);
        time_ch.resize(number_of_nodes);
        hf_dis.resize(number_of_nodes);
        CH_order.resize(number_of_nodes);
        dis.resize(number_of_nodes);
        pos.resize(number_of_nodes);
        X.resize(number_of_nodes);
        fi_dis.resize(number_of_nodes);
        fa.resize(number_of_nodes);
        ancestor.resize(number_of_nodes);
        dep.resize(number_of_nodes);
        H2H_order.resize(number_of_nodes);
        for (int i = 0; i < number_of_nodes; i++) {
            file >> ix >> dx >> dy;
            lon[i] = dx;
            lat[i] = dy;
            nodes[i].osmid=i;
        }
        file.close();
    } else {
        bo=0;
        std::cout << "Unable to open file." << std::endl;
    }
    if(!bo)return;
    std::cout<<"ready nodes<<'\n";

    //get the information of streets
    file.open(graph_file);

    if (file.is_open()) {
        file >> number_of_line;
        for (int i = 0; i < number_of_line; i++) {
            file>>ix>>iy>>dx;
            if(exist_line[ix][iy])
            {
                for(auto& edge:graph[ix])
                {
                    if(edge.first==iy)
                    {
                        edge.second=std::min(edge.second,(int)dx);
                        break;
                    }
                }
                for(auto& edge:graph[iy])
                {
                    if(edge.first==ix)
                    {
                        edge.second=std::min(edge.second,(int)dx);
                        break;
                    }
                }
                continue;
            }
            else
            {
                graph[ix].push_back(std::make_pair(iy,dx));
                graph[iy].push_back(std::make_pair(ix,dx));
                nodes[ix].degree++;
                nodes[iy].degree++;
            }
            
            exist_line[ix][iy]=exist_line[iy][ix]=1;
        }
        file.close();
    } else {
        bo=0;
        std::cout << "Unable to open file." << std::endl;
    }
    if(!bo)return;
    std::cout<<"ready graph<<'\n";
    //readin a social graph
    // std::ifstream file(graph_file);
    // int ix,iy,dx;
    // if (file.is_open()) {
    //     file >> number_of_line;
    //     lon.resize(number_of_nodes);
    //     lat.resize(number_of_nodes);
    //     pre.resize(number_of_nodes);
    //     distance.resize(number_of_nodes);
    //     graph.resize(number_of_nodes);
    //     predict.resize(number_of_nodes);
    //     frontdis.resize(number_of_nodes);
    //     backdis.resize(number_of_nodes);
    //     nodes.resize(number_of_nodes);
    //     Labels.resize(number_of_nodes);
    //     order.resize(number_of_nodes);
    //     contracted_nodes.resize(number_of_nodes);
    //     time_ch.resize(number_of_nodes);
    //     hf_dis.resize(number_of_nodes);
    //     CH_order.resize(number_of_nodes);
    //     dis.resize(number_of_nodes);
    //     pos.resize(number_of_nodes);
    //     X.resize(number_of_nodes);
    //     fi_dis.resize(number_of_nodes);
    //     fa.resize(number_of_nodes);
    //     ancestor.resize(number_of_nodes);
    //     dep.resize(number_of_nodes);
    //     H2H_order.resize(number_of_nodes);
    //     for (int i = 0; i < number_of_nodes; i++) nodes[i].osmid=i;
    //     for (int i = 0; i < number_of_line; i++) {
    //         file>>ix>>iy>>dx;
    //         graph[ix].push_back(std::make_pair(iy,dx));
    //         graph[iy].push_back(std::make_pair(ix,dx));
    //         nodes[ix].degree++;
    //         nodes[iy].degree++;
    //     }
    //     file.close();
    // } else {
    //     bo=0;
    //     std::cout << "Unable to open file." << std::endl;
    // }
    // if(!bo)return;
    // std::cout<<"!";
    preprocess2Hop();
    // work();
    // self_test();
}

void Graph::work()
{
    cnt=0;
    std::ifstream infile("Facebook.query");
    infile>>number_of_query;
    test.resize(number_of_query);
    std::cout<<number_of_query<<std::endl;
    std::ofstream outfile("Facebook.answerPLL");
    clock_t start = clock();
    for(int i=0;i<number_of_query;i++)
    {
        int s,t;
        infile>>s>>t;
        outfile<<q2Hop(s,t)<<std::endl;
        // if(i%100==99)
        // {
        //     clock_t end = clock();
        //     double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
        //     std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        //     start=end;
        // }
    }
    infile.close();
}

void Graph::preprocess2Hop()
{
    clock_t start=clock();
    sort(nodes.begin(),nodes.end(),comparenodes);
    for(int i=0;i < number_of_nodes; i++)
    {
        order[nodes[i].osmid]=i;
    }
    std::cout<<"ready order"<<'\n';
    for (int i = 0; i < number_of_nodes; ++i) {
        prunedBFS(nodes[i].osmid);
        std::cout<<"BFS once"<<'\n';
        // if(i%(number_of_nodes/1000)==0)
        // {
        //     std::cout<<nodes[i].degree<<" "<<nodes[i].osmid<<std::endl;
        //     clock_t end = clock();
        //     double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
        //     std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        // }
    }
    std::cout<<"ready BFS"<<'\n';
    for(int i=0;i<number_of_nodes;i++)
    {
        std::sort(Labels[i].begin(),Labels[i].end(),compareLabel);
    }
}
void Graph::prunedBFS(int start) {
    heap q(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
    }
    distance[start]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    q.update(start,0);
    while(q.size())
    {
        //cnt++;
        int element,key;
        q.extract_min(element,key);
        if(vis[element])continue;
        vis[element]=1;
        Labels[start].push_back({element, distance[element]});
        Labels[element].push_back({start, distance[element]});
        for(auto edge : graph[element])
        {
            int v=edge.first;
            if(distance[v]>distance[element]+edge.second)
            {
                distance[v]=distance[element]+edge.second;
                //pre[v]=element;
                bool prune= false;
                for (const auto& p : Labels[v]) {
                    int mid = p.vertex;
                    int distMidToV = p.distance;
                    if (distance[start] + distance[mid] + distMidToV <= distance[v]) {
                        prune = true;
                        break;
                    }
                }
                if (!prune) {
                    q.update(v,distance[v]);
                }
            }
        }
    }
    
}

int Graph::q2Hop(int u, int v) {
    int min_distance = std::numeric_limits<int>::max();
    const auto& labels_u = Labels[u];
    const auto& labels_v = Labels[v];
    size_t i = 0, j = 0;
    while (i < labels_u.size() && j < labels_v.size()) {
        if (labels_u[i].vertex == labels_v[j].vertex) {
            int dis=labels_u[i].distance + labels_v[j].distance;
            if(min_distance>dis)
            {
                min_distance=dis;
                //std::cout<<labels_u[i].vertex<<std::endl;
            }
            ++i;
            ++j;
        } else if (labels_u[i].vertex < labels_v[j].vertex) {
            ++i;
        } else {
            ++j;
        }
    }
    return min_distance;
}
std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::PLL_way(int s, int t, int dis) {
    std::vector<std::pair<double, double>> path;
    std::vector<std::pair<double, double>> search;

    while (s != t) {
        path.push_back({lat[s],lon[s]});
        search.push_back({lat[s],lon[s]});
        for (auto edge : graph[s]) {
            int v = edge.first, w = edge.second;
            if (w + q2Hop(v, t) == dis) {
                s = v;
                dis -= w;
                break;
            }
        }
    }

    path.push_back({lat[s],lon[s]});
    search.push_back({lat[s],lon[s]});

    return {path, search};
}