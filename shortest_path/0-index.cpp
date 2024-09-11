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
    readin(node_file, edge_file);
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
    if (!bo) return;
}


// void Graph::work()
// {
//     cnt=0;
//     std::ifstream infile("GZ.query");
//     infile>>number_of_query;
//     test.resize(number_of_query);
//     std::cout<<number_of_query<<std::endl;
//     std::ofstream outfile("GZ.answerAstar");
//     clock_t start = clock();
//     for(int i=0;i<number_of_query;i++)
//     {
//         int s,t;
//         infile>>s>>t;
//         outfile<<Astar(s,t)<<std::endl;
//         if(i%100==99)
//         {
//             clock_t end = clock();
//             double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
//             std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
//             start=end;
//         }
//     }
//     infile.close();
// }


int Graph::Astar(int s,int t)
{
    //std::cout<<s<<" "<<t<<std::endl;
    //pretreatment
    heap q(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
        predict[i]=0;
    }
    distance[s]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    //using dijistra and A*
    q.update(s,0);
    while(q.size())
    {
        //cnt++;
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
                pre[v]=element;
                //std::cout<<v<<std::endl;
                q.update(v,distance[v]+real_dis(lat[v],lon[v],lat[t],lon[t]));
                //q.update(v,distance[v]);
                //std::cout<<v<<std::endl;
            }
        }
    }
    //std::cout<<t<<std::endl;
    return distance[t];
}

int Graph::dijkstra(int s,int t)
{
    //pretreatment
    heap q(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
    }
    distance[s]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    //using dijistra and A*
    q.update(s,0);
    while(q.size())
    {
        int element,key;
        q.extract_min(element,key);
        if(element==t)break;
        vis[element]=1;
        for(auto edge : graph[element])
        {
            int v=edge.first;
            if(vis[v])continue;
            if(distance[v]>distance[element]+edge.second)
            {
                distance[v]=distance[element]+edge.second;
                pre[v]=element;
                q.update(v,distance[v]);
                vis[v]=0;
            }
        }
    }
    if(t==number_of_nodes)return 0;
    return distance[t];
}
// std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::dijkstra_way(int s, int t) {
//     heap q(number_of_nodes);
//     for (int i = 0; i < number_of_nodes; i++) {
//         distance[i] = 1e9;
//     }
//     distance[s] = 0;
//     std::vector<int> vis;
//     vis.resize(number_of_nodes, 0);

//     q.update(s, 0);
//     while (q.size()) {
//         int element, key;
//         q.extract_min(element, key);
//         if (element == t) break;
//         vis[element] = 1;
//         for (auto edge : graph[element]) {
//             int v = edge.first;
//             if (vis[v]) continue;
//             if (distance[v] > distance[element] + edge.second) {
//                 distance[v] = distance[element] + edge.second;
//                 pre[v] = element;
//                 q.update(v, distance[v]);
//                 vis[v] = 0;
//             }
//         }
//     }

//     std::vector<std::pair<double, double>> path;
//     std::vector<std::pair<double, double>> search;
//     if(distance[t]==1e9){
//         return(std::make_pair(path,search));
//     }else{
//         for (int at = t;at!=-1;at=pre[at]){
//             path.push_back(std::make_pair(lat[at],lon[at]));
//         }
//         std::reverse(path.begin(),path.end());
//         for(int i = 0;i < number_of_nodes;i++){
//             if (vis[i]){
//                 search.push_back(std::make_pair(lat[i],lon[i]));
//             }
//         }
//     }
//     return(std::make_pair(path,search));
// }

std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::dijkstra_way(int s, int t) {
    // std::cout << "Starting Dijkstra way from " << s << " to " << t << std::endl;

    heap q(number_of_nodes);
    for (int i = 0; i < number_of_nodes; i++) {
        distance[i] = 1e9;
    }
    distance[s] = 0;
    std::vector<int> vis;
    vis.resize(number_of_nodes, 0);

    q.update(s, 0);
    while (q.size()) {
        int element, key;
        q.extract_min(element, key);
        // std::cout << "Extracted node " << element << " with key " << key << std::endl;

        if (element == t) break;
        vis[element] = 1;

        for (auto edge : graph[element]) {
            int v = edge.first;
            if (vis[v]) continue;

            if (distance[v] > distance[element] + edge.second) {
                distance[v] = distance[element] + edge.second;
                pre[v] = element;
                q.update(v, distance[v]);
                vis[v] = 0;
                // std::cout << "Updated distance for node " << v << " to " << distance[v] << std::endl;
            }
        }
    }

    std::vector<std::pair<double, double>> path;
    std::vector<std::pair<double, double>> search;

    if (distance[t] == 1e9) {
        std::cerr << "No valid path found from " << s << " to " << t << std::endl;
        return std::make_pair(path, search);
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

    // std::cout << "Dijkstra way found a path of length: " << path.size() << std::endl;

    return std::make_pair(path, search);
}


int Graph::Bidijkstra(int s,int t)
{
    heap fronth(number_of_nodes),backh(number_of_nodes);
    std::vector<bool>closed_s,closed_t;
    closed_s.resize(number_of_nodes,false);
    closed_t.resize(number_of_nodes,false);
    int answer=1e9;
    for(int i=0;i<number_of_nodes;i++)
    {
        frontdis[i]=backdis[i]=1e9;
    }
    frontdis[s]=backdis[t]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    fronth.update(s,0);
    backh.update(t,0);
    while(fronth.size()&&backh.size())
    {
        cnt+=2;
        int element_s,key_s,element_t,key_t;
        fronth.extract_min(element_s,key_s);
        backh.extract_min(element_t,key_t);
        closed_s[element_s]=true;
        if(closed_t[element_s])break;
        for(auto edge : graph[element_s])
        {
            int v=edge.first;
            //std::cout<<v<<std::endl;
            if(closed_t[v])answer=std::min(answer,frontdis[element_s]+edge.second+backdis[v]);
            if(!closed_s[v]&& frontdis[v]>frontdis[element_s]+edge.second)
            {
                frontdis[v]=frontdis[element_s]+edge.second;
                fronth.update(v,frontdis[v]);
                //std::cout<<v<<std::endl;
            }
        }

        closed_t[element_t]=true;
        if(closed_s[element_t])break;
        for(auto edge : graph[element_t])
        {
            int v=edge.first;
            if(closed_s[v])answer=std::min(answer,frontdis[v]+edge.second+backdis[element_t]);
            if(!closed_t[v]&& backdis[v]>backdis[element_t]+edge.second)
            {
                backdis[v]=backdis[element_t]+edge.second;
                backh.update(v,backdis[v]);
            }
        }
    }
    //std::cout<<t<<std::endl;
    return answer;

}

std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::bi_dijkstra_way(int s,int t){
    heap fronth(number_of_nodes), backh(number_of_nodes);
    std::vector<bool> closed_s(number_of_nodes, false), closed_t(number_of_nodes, false);
    int answer = 1e9;
    for (int i = 0; i < number_of_nodes; i++) {
        frontdis[i] = backdis[i] = 1e9;
    }
    // std::cout<<"preokey"<<'\n';
    frontdis[s] = backdis[t] = 0;
    std::vector<int> pre_s(number_of_nodes, -1), pre_t(number_of_nodes, -1);
    fronth.update(s, 0);
    backh.update(t, 0);
    bool meet = false;
    int meet_node = -1;
    
    while (fronth.size() && backh.size()) {
        int element_s, key_s, element_t, key_t;
        fronth.extract_min(element_s, key_s);
        backh.extract_min(element_t, key_t);

        closed_s[element_s] = true;
        closed_t[element_t] = true;
        if (closed_t[element_s]) {
            meet = true;
            meet_node = element_s;
            break;
        }
        for (auto edge : graph[element_s]) {
            int v = edge.first;
            if (closed_t[v]) answer = std::min(answer, frontdis[element_s] + edge.second + backdis[v]);
            if (!closed_s[v] && frontdis[v] > frontdis[element_s] + edge.second) {
                frontdis[v] = frontdis[element_s] + edge.second;
                pre_s[v] = element_s;
                fronth.update(v, frontdis[v]);
            }
        }

        if (closed_s[element_t] && closed_t[element_t]) {
            meet = true;
            meet_node = element_t;
            // std::cout<<meet_node;
            break;
            
        }
        
        if (closed_s[element_t]) {
            meet = true;
            meet_node = element_s;
            break;
        }
        for (auto edge : graph[element_t]) {
            int v = edge.first;
            if (closed_s[v]) answer = std::min(answer, frontdis[v] + edge.second + backdis[element_t]);
            if (!closed_t[v] && backdis[v] > backdis[element_t] + edge.second) {
                backdis[v] = backdis[element_t] + edge.second;
                pre_t[v] = element_t;
                backh.update(v, backdis[v]);
            }
        }
    }
    // std::cout<<"disok"<<"ans:"<<answer<<std::endl;
    std::vector<std::pair<double, double>> path;
    std::vector<std::pair<double, double>> search;
    
    if (!meet) {
        return std::make_pair(path, search);
    } else {
       
        for (int at = meet_node; at != -1; at = pre_s[at]) {
            path.push_back(std::make_pair(lat[at], lon[at]));
        }
        // std::cout<<"wayok";
        std::reverse(path.begin(), path.end());

        for (int at = pre_t[meet_node]; at !=-1; at = pre_t[at]) {
            path.push_back(std::make_pair(lat[at], lon[at]));
        }
        for (int i = 0; i < number_of_nodes;i++) {
            if (closed_s[i] || closed_t[i]) {
                search.push_back(std::make_pair(lat[i], lon[i]));
            }
        }

        return std::make_pair(path, search);
    }
}
std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::Astar_way(int s,int t){
    heap q(number_of_nodes);
    for(int i=0;i<number_of_nodes;i++)
    {
        distance[i]=1e9;
        predict[i]=0;
    }
    distance[s]=0;
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
    //using dijistra and A*
    q.update(s,0);
    while(q.size())
    {
        //cnt++;
        int element,key;
        q.extract_min(element,key);
        if(element==t)break;
        vis[element]=1;
        //std::cout<<element<<std::endl;
        for(auto edge : graph[element])
        {
            int v=edge.first;
            //std::cout<<v<<std::endl;
            if(distance[v]>distance[element]+edge.second)
            {
                distance[v]=distance[element]+edge.second;
                pre[v]=element;
                //std::cout<<v<<std::endl;
                q.update(v,distance[v]+real_dis(lat[v],lon[v],lat[t],lon[t]));
                //q.update(v,distance[v]);
                //std::cout<<v<<std::endl;
            }
        }
    }
    std::vector<std::pair<double, double>> path;
    std::vector<std::pair<double, double>> search;
    if(distance[t]==1e9){
        return(std::make_pair(path,search));
    }else{
        for (int at = t;at!=-1;at=pre[at]){
            path.push_back(std::make_pair(lat[at],lon[at]));
        }
        std::reverse(path.begin(),path.end());
        for(int i = 0;i < number_of_nodes;i++){
            if (vis[i]){
                search.push_back(std::make_pair(lat[i],lon[i]));
            }
        }
    }
    return(std::make_pair(path,search));
    
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
