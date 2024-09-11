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

Graph::Graph(const std::string& node_file, const std::string& edge_file) 
{
    readin(node_file,edge_file);
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

void Graph::readin(const std::string& node_file, const std::string& edge_file)
{
    int bo = 1;
    exist_line.clear();
    // Read in a traffic graph
    std::ifstream file(node_file);
    long long ix, iy;
    double dx, dy;

    // Get the information of nodes
    if (file.is_open()) {
        file >> number_of_nodes;
        lon.resize(number_of_nodes);
        lat.resize(number_of_nodes);
        pre.resize(number_of_nodes,-1);
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
            lon[i] = dx ;
            lat[i] = dy ;
            nodes[i].osmid = i;
            // std::cout<<dx<<dy;
            
        }
        file.close();
    } else {
        bo = 0;
        std::cout << "Unable to open node file." << std::endl;
    }
    if (!bo) return;

    // Get the information of streets
    file.open(edge_file);

    if (file.is_open()) {
        file >> number_of_line;
        for (int i = 0; i < number_of_line; i++) {
            file >> ix >> iy >> dx;
            // std::cout<<ix<<iy<<std::endl;
            if (exist_line[ix][iy]) {
                for (auto& edge : graph[ix]) {
                    if (edge.first == iy) {
                        edge.second = std::min(edge.second, (int)dx);
                        break;
                    }
                }
                for (auto& edge : graph[iy]) {
                    if (edge.first == ix) {
                        edge.second = std::min(edge.second, (int)dx);
                        break;
                    }
                }
                continue;
            } else {
                graph[ix].push_back(std::make_pair(iy, dx));
                graph[iy].push_back(std::make_pair(ix, dx));
                // std::cout<<ix<<" "<<iy<<'/n';
                nodes[ix].degree++;
                nodes[iy].degree++;
            }

            exist_line[ix][iy] = exist_line[iy][ix] = 1;
        }
        file.close();
    } else {
        bo = 0;
        std::cout << "Unable to open edge file." << std::endl;
    }
    std::cout<<"!";
    Preprocessed_path.resize(20);
    for(int i=0;i<20;i++)
    {
        preprocess(i,number_of_nodes/20*i);
    }
    if (!bo) return;
}

    //readin a social graph
    /*std::ifstream file("Facebook");
    int ix,iy,dx;
    if (file.is_open()) {
        file >> number_of_nodes>> number_of_line;
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
        for (int i = 0; i < number_of_nodes; i++) nodes[i].osmid=i;
        for (int i = 0; i < number_of_line; i++) {
            file>>ix>>iy>>dx;
            graph[ix].push_back(std::make_pair(iy,dx));
            graph[iy].push_back(std::make_pair(ix,dx));
            nodes[ix].degree++;
            nodes[iy].degree++;
        }
        file.close();
    } else {
        bo=0;
        std::cout << "Unable to open file." << std::endl;
    }
    if(!bo)return;*/
    
    // work();


void Graph::work()
{
    cnt=0;
    std::ifstream infile("GZ.query");
    infile>>number_of_query;
    test.resize(number_of_query);
    std::cout<<number_of_query<<std::endl;
    std::ofstream outfile("GZ.answerALT");
    clock_t start = clock();
    for(int i=0;i<number_of_query;i++)
    {
        int s,t;
        infile>>s>>t;
        outfile<<ALT(s,t)<<std::endl;
        if(i%100==99)
        {
            clock_t end = clock();
            double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
            std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
            start=end;
        }
    }
    infile.close();
}

const double EARTH_RADIUS_KM = 6371.0;
int Graph::real_dis(double lat1, double lon1, double lat2, double lon2)
{
    double lat1Rad = degreesToRadians(lat1);
    double lon1Rad = degreesToRadians(lon1);
    double lat2Rad = degreesToRadians(lat2);
    double lon2Rad = degreesToRadians(lon2);

    double dLat = lat2Rad - lat1Rad;
    double dLon = lon2Rad - lon1Rad;

    double a = std::sin(dLat / 2) * std::sin(dLat / 2) +
               std::cos(lat1Rad) * std::cos(lat2Rad) *
               std::sin(dLon / 2) * std::sin(dLon / 2);

    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    return (int)(EARTH_RADIUS_KM * c*1000);
}

int Graph::query_h(int node,int t)
{
    //return 0;
    if(node==t)return 0;
    if(predict[node])return predict[node];
    int mx=real_dis(lat[node], lon[node], lat[t], lon[t]);
    //std::cout<<"h"<<std::endl;
    for(int i=0;i<20;i++)mx=std::max(mx,Preprocessed_path[i][node]-Preprocessed_path[i][t]);
    return mx;
}

void Graph::preprocess(int order,int s)
{
    //std::cout<<order<<" "<<s<<std::endl;
    std::priority_queue<std::pair<int, int> >q;
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
    }
    distance[s]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    //using dijistra
    q.push(std::make_pair(0,s));
    while(q.size())
    {
        int x=q.top().second;q.pop();
        if(vis[x])continue;
        vis[x]=1;
        //std::cout<<x<<std::endl;
        for(auto edge : graph[x])
        {
            int v=edge.first;
            if(distance[v]>distance[x]+edge.second)
            {
                distance[v]=distance[x]+edge.second;
                // pre[v]=x;
                q.push(std::make_pair(-distance[v],v));
            }
        }
    }
    Preprocessed_path[order].resize(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        Preprocessed_path[order][i]=distance[i];
    }
}

int Graph::ALT(int s,int t)
{
    //pretreatment
    heap q(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
        predict[i]=0;
    }
    distance[s]=0;
    // std::vector<int>vis;
    // vis.resize(number_of_nodes,0);
    //using dijistra and A*
    q.update(s,0);
    while(q.size())
    {
        int element,key;
        q.extract_min(element,key);
        if(element==t)break;
        //std::cout<<element<<std::endl;
        for(auto edge : graph[element])
        {
            int v=edge.first;
            //std::cout<<v<<std::endl;
            if(distance[v]>distance[element]+edge.second)
            {
                distance[v]=distance[element]+edge.second;
            
                //std::cout<<v<<std::endl;
                q.update(v,distance[v]+query_h(v,t));
                //q.update(v,distance[v]);
                //std::cout<<v<<std::endl;
            }
        }
    }
    //std::cout<<t<<std::endl;
    return distance[t];
}
std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::ALT_way(int s, int t) {
    // Pretreatment
    heap q(number_of_nodes);
    for (int i = 0; i < number_of_nodes; i++) {
        distance[i] = 1e9;
        predict[i] = 0;
    }
    distance[s] = 0;
    std::vector<int> vis(number_of_nodes, 0);
    std::cout << "why" << std::endl;

    // Using Dijkstra and A*
    q.update(s, 0);
    while (q.size()) {
        int element, key;
        q.extract_min(element, key);
        

        if (element == t) {

            break;
        }
        vis[element] = 1;

        for (auto edge : graph[element]) {
            int v = edge.first;
            double new_dist = distance[element] + edge.second;
            // std::cout << "Checking edge to: " << v << " with new distance: " << new_dist << std::endl;

            if (distance[v] > new_dist) {
                distance[v] = new_dist;
                pre[v] = element;
                double heuristic = query_h(v, t);
                // std::cout << "Updating node: " << v << " with distance: " << distance[v] << " and heuristic: " << heuristic << std::endl;
                q.update(v, distance[v] + heuristic);
            }
        }
    }

    // Prepare path and search vectors
    std::vector<std::pair<double, double>> path;
    std::vector<std::pair<double, double>> search;

    if (distance[t] == 1e9) {
        std::cout << "No valid path found." << std::endl;
        return std::make_pair(path, search); // No valid path found
    } else {
        for (int at = t; at != -1; at = pre[at]) {
            path.push_back(std::make_pair(lat[at], lon[at]));
        }
        std::reverse(path.begin(), path.end());

        for (int i = 0; i < number_of_nodes; i++) {
            if (vis[i]) {
                search.push_back(std::make_pair(lat[i], lon[i]));
            }
        }
    }
    return std::make_pair(path, search);
}