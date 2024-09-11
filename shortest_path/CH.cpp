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

void Graph::readin(const std::string&Filename_N,const std::string &Filename_G) 
{
    int bo=1;
    exist_line.clear();
    //readin a traffic graph
    std::ifstream file(Filename_N);
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
            // lon[i] = dy;
            // lat[i] = dx;
            nodes[i].osmid=i;
        }
        file.close();
    } else {
        bo=0;
        std::cout << "Unable to open file." << std::endl;
    }
    if(!bo)return;

    //get the information of streets
    file.open(Filename_G);

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
    // if (load_preprocessed_data(CH_order, shortcut_map)) {
    //     // 成功加载预处理数据，跳过预处理步骤
    //     std::cout << "Preprocessed data loaded successfully, skipping preprocessing." << std::endl;
    //     return;
    // } else {
    //     std::cout << "No preprocessed data found, starting preprocessing." << std::endl;
    //     CH_preprocess();
    // }
    CH_preprocess();
   
    // CH_preprocess();
    // work();
}

// void Graph::work()
// {
//     cnt=0;
//     std::ifstream infile("GZ.query");
//     infile>>number_of_query;
//     test.resize(number_of_query);
//     std::cout<<number_of_query<<std::endl;
//     std::ofstream outfile("GZ.answerCH");
//     clock_t start = clock();
//     for(int i=0;i<number_of_query;i++)
//     {
//         int s,t;
//         infile>>s>>t;
//         outfile<<CH(s,t)<<std::endl;
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

void Graph::CH_preprocess()
{
    heap q(number_of_nodes);
    for (int i = 0; i < number_of_nodes; i++) {
        q.update(i, nodes[i].degree);
    }
    // clock_t start=clock();
    while (!q.empty()) {
        int node, priority;
        q.extract_min(node, priority);
        CH_order[node]=cnt;
        cnt++;
        contracted_nodes[node]=1;
        if(nodes[node].degree==1)
        {
            for (const auto& neighbora : graph[node])
            {
                int v=neighbora.first;
                if(contracted_nodes[v])continue;
                nodes[v].degree--;
                q.update(v,nodes[v].degree);
            }
            continue;
        }
        for (const auto& neighbora : graph[node]) {
            int v = neighbora.first;
            if (contracted_nodes[v]) continue;

            for (const auto& neighborb : graph[node]) {
                int w = neighborb.first;
                int shortcut = neighbora.second + neighborb.second;
                if (contracted_nodes[w]) continue;
                if (v == w) continue;

                if (exist_line[v][w]) {
                    bool updated = false;
                    for (auto& edge : graph[v]) {
                        if (edge.first == w) {
                            if (shortcut < edge.second) {
                                edge.second = shortcut;
                                updated = true;
                            }
                            break;
                        }
                    }
                    for (auto& edge : graph[w]) {
                        if (edge.first == v) {
                            if (shortcut < edge.second) {
                                edge.second = shortcut;
                                updated = true;
                            }
                            break;
                        }
                    }

                    if (updated) {
                        std::vector<int> part;
                        if (shortcut_map.find({v, node}) != shortcut_map.end()) {
                            part = shortcut_map[{v, node}];
                        } else {
                            part.push_back(v);
                        }

                        part.push_back(node);
                        if (shortcut_map.find({node, w}) != shortcut_map.end()) {
                            std::vector<int> subpart = shortcut_map[{node, w}];
                            part.insert(part.end(), subpart.begin(), subpart.end());
                        } else {
                            part.push_back(w);
                        }

                        shortcut_map[{v, w}] = part;
                        shortcut_map[{w, v}] = std::vector<int>(part.rbegin(), part.rend());
                    }

                } else {
        
                    graph[v].push_back(std::make_pair(w, shortcut));
                    graph[w].push_back(std::make_pair(v, shortcut));
                    nodes[v].degree++;
                    nodes[w].degree++;
                    exist_line[v][w] = exist_line[w][v] = 1;

                    std::vector<int> part;
                    if (shortcut_map.find({v, node}) != shortcut_map.end()) {
                        part = shortcut_map[{v, node}];
                    } else {
                        part.push_back(v);
                    }

                    part.push_back(node);
                    if (shortcut_map.find({node, w}) != shortcut_map.end()) {
                        std::vector<int> subpart = shortcut_map[{node, w}];
                        part.insert(part.end(), subpart.begin(), subpart.end());
                    } else {
                        part.push_back(w);
                    }

                    shortcut_map[{v, w}] = part;
                    shortcut_map[{w, v}] = std::vector<int>(part.rbegin(), part.rend());
                }
            }
        }
        
        for (const auto& neighbora : graph[node])
        {
            int v=neighbora.first;
            if(contracted_nodes[v])continue;
            nodes[v].degree--;
            q.update(v,nodes[v].degree);
        }
        
        // if (cnt % (number_of_nodes / 100) == 0) {
        //     std::cout << priority << " " << node <<" "<<cnt<< std::endl;
        //     clock_t end = clock();
        //     double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
        //     std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        // }
    }
    // std::cout<<"done preprocess"<<std::endl;
    for(int i=0;i<number_of_nodes;i++)contracted_nodes[i]=0;
    // std::cout<<"process ready"<<'\n';
    
    // save_preprocessed_data(CH_order, shortcut_map);

}





int Graph::CH(int s,int t)
{
    if(s==t)return 0;
    //std::cout<<s<<" "<<t<<std::endl;
    heap fronth(number_of_nodes),backh(number_of_nodes);
    std::vector<bool>closed_s,closed_t;
    closed_s.resize(number_of_nodes,false);
    closed_t.resize(number_of_nodes,false);

    std::vector<int> pre_s,pre_t;
    pre_s.resize(number_of_nodes,-1);
    pre_t.resize(number_of_nodes,-1);
    std::vector<bool> stc_s,stc_t;
    stc_s.resize(number_of_nodes,false);
    stc_t.resize(number_of_nodes,false);
    // std::cout<<"prepare ready"<<'\n';
    // std::vector<int> pre_s;
    // std::vector<int>pre_t;
    int answer=1e9;
    for(int i=0;i<number_of_nodes;i++)
    {
        frontdis[i]=backdis[i]=1e9;
    }
    frontdis[s]=backdis[t]=0;
    //using dijistra and A*
    fronth.update(s,0);
    backh.update(t,0);
    int ti=0;
    while(fronth.size()||backh.size())
    {
        // if(fronth.top()+backh.top()>=answer)break;
        //if(ti>=10)break;
        
        cnt+=2;
        int element_s,key_s,element_t,key_t;

        if(fronth.size())
        {   
            fronth.extract_min(element_s,key_s);
            closed_s[element_s]=true;
            answer=std::min(answer,frontdis[element_s]+backdis[element_s]);
            //std::cout<<element<<std::endl;
            for(auto edge : graph[element_s])
            {   
                int v=edge.first;
                if(CH_order[v]<CH_order[element_s])continue;
                //std::cout<<v<<std::endl;
                if(closed_t[v]) answer=std::min(answer,frontdis[element_s]+edge.second+backdis[v]);
                if(!closed_s[v]&& frontdis[v]>frontdis[element_s]+edge.second)
                {   
                    frontdis[v]=frontdis[element_s]+edge.second;
                    // pre_s[v]=element_s;
                    // stc_s[v]=shortcut_map.find({element_s,v})!=shortcut_map.end();
                    fronth.update(v,frontdis[v]);
                    //std::cout<<v<<std::endl;
                }
            }
            if(closed_t[element_s])ti++;
        }
        if(backh.size())
        {      
            
            backh.extract_min(element_t,key_t);
            closed_t[element_t]=true;
            answer=std::min(answer,frontdis[element_t]+backdis[element_t]);
            for(auto edge : graph[element_t])
            {   
                int v=edge.first;
                if(CH_order[v]<CH_order[element_t])continue;
                //std::cout<<v<<std::endl;
                if(closed_s[v])answer=std::min(answer,frontdis[v]+edge.second+backdis[element_t]);
                if(!closed_t[v]&& backdis[v]>backdis[element_t]+edge.second)
                {  
                    backdis[v]=backdis[element_t]+edge.second;
                    // pre_t[v]=element_t;
                    // stc_t[v]=shortcut_map.find({element_t,v})!=shortcut_map.end();
                    backh.update(v,backdis[v]);
                    //std::cout<<v<<std::endl;
                }
            }
            if(closed_s[element_t])ti++;
        }
        
    }
    //std::cout<<t<<std::endl;
    return answer;
}
std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::CH_way(int s, int t) {
    heap fronth(number_of_nodes), backh(number_of_nodes);
    std::vector<bool> closed_s(number_of_nodes, false), closed_t(number_of_nodes, false);
    std::vector<int>  pre_s(number_of_nodes, -1), pre_t(number_of_nodes, -1);
    std::vector<bool> stc_s(number_of_nodes, false), stc_t(number_of_nodes, false);

    int answer = 1e9;
    int meet = -1; 
    for (int i = 0; i < number_of_nodes; i++) {
        frontdis[i] = backdis[i] = 1e9;
    }
    frontdis[s] = backdis[t] = 0;

    fronth.update(s, 0);
    backh.update(t, 0);
    int ti = 0;

    std::vector<std::pair<double, double>> search_points; 

    while (fronth.size() || backh.size()) {
        cnt += 2;
        int element_s = -1, key_s = -1, element_t = -1, key_t = -1;

        if (fronth.size()) {
            fronth.extract_min(element_s, key_s);
            closed_s[element_s] = true;
            answer = std::min(answer, frontdis[element_s] + backdis[element_s]);
            search_points.push_back({lat[element_s], lon[element_s]}); 
            for (auto edge : graph[element_s]) {
                int v = edge.first;
                if (CH_order[v] < CH_order[element_s]) continue;

                if (closed_t[v] && frontdis[element_s] + edge.second + backdis[v] < answer) {
                    answer = std::min(answer, frontdis[element_s] + edge.second + backdis[v]);
                    meet = v; 
                }
                if (!closed_s[v] && frontdis[v] > frontdis[element_s] + edge.second) {
                    frontdis[v] = frontdis[element_s] + edge.second;
                    pre_s[v] = element_s;
                    stc_s[v] = shortcut_map.find({element_s, v}) != shortcut_map.end();
                    fronth.update(v, frontdis[v]);
                }
            }

            if (closed_t[element_s]) ti++;
        }

        if (backh.size()) {
            backh.extract_min(element_t, key_t);
            closed_t[element_t] = true;
            answer = std::min(answer, frontdis[element_t] + backdis[element_t]);

            search_points.push_back({lat[element_t], lon[element_t]}); 


            for (auto edge : graph[element_t]) {
                int v = edge.first;
                if (CH_order[v] < CH_order[element_t]) continue;
                if (closed_s[v]&& frontdis[v] + edge.second + backdis[element_t] < answer) {
                    answer = std::min(answer, frontdis[v] + edge.second + backdis[element_t]);
                    meet = v;
                }
                if (!closed_t[v] && backdis[v] > backdis[element_t] + edge.second) {
                    backdis[v] = backdis[element_t] + edge.second;
                    pre_t[v] = element_t;
                    stc_t[v] = shortcut_map.find({element_t, v}) != shortcut_map.end();
                    backh.update(v, backdis[v]);
                }
            }
            if (closed_s[element_t]) ti++;
        }
    }

    if (meet == -1) {
        std::cout << "No meeting point found!" << std::endl;
        return {{}, search_points}; 
    }

    std::vector<std::pair<double, double>> path;

    // // 调试输出 - 路径回溯开始
    // std::cout << "Tracing path from s to meet, meet = " << meet << std::endl;

    int current = meet;

    // 回溯从 meet 到 s 的路径
    while (current != s) {
        // std::cout << "Current node in path (s to meet): " << current << std::endl;
        if (stc_s[current]) {
            std::vector<int>& inside = shortcut_map[{pre_s[current], current}];
            // std::cout << "Expanding shortcut from " << pre_s[current] << " to " << current << std::endl;
            for (auto it = inside.rbegin(); it != inside.rend(); ++it) {
                path.push_back({lat[*it], lon[*it]});
            }
        } else {
            path.push_back({lat[current], lon[current]});
        }
        current = pre_s[current];
    }
    path.push_back({lat[s],lon[s]});

    std::reverse(path.begin(), path.end());

    // std::cout << "Tracing path from meet to t, t = " << t << std::endl;

    current = meet;
    while (current != t) {
        // std::cout << "Current node in path (meet to t): " << current << std::endl;
        if (stc_t[current]) {
            std::vector<int>& inside = shortcut_map[{current, pre_t[current]}];
            // std::cout << "Expanding shortcut from " << current << " to " << pre_t[current] << std::endl;
            for (auto it = inside.begin(); it != inside.end(); ++it) {
                path.push_back({lat[*it], lon[*it]});
            }
        } else {
            path.push_back({lat[pre_t[current]], lon[pre_t[current]]});
        }
        current = pre_t[current];
    }

    path.push_back({lat[t], lon[t]});


    // std::cout << "Finished path tracing." << std::endl;

    return {path, search_points};
}

// std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> Graph::CH_way(int s, int t) {
//     heap fronth(number_of_nodes), backh(number_of_nodes);
//     std::vector<bool> closed_s(number_of_nodes, false);
//     std::vector<bool> closed_t(number_of_nodes, false);

//     std::vector<int> pre_s(number_of_nodes, -1);
//     std::vector<int> pre_t(number_of_nodes, -1);
//     std::vector<bool> stc_s(number_of_nodes, false);
//     std::vector<bool> stc_t(number_of_nodes, false);

//     int answer = 1e9;
//     for (int i = 0; i < number_of_nodes; i++) {
//         frontdis[i] = backdis[i] = 1e9;
//     }
//     frontdis[s] = 0;
//     backdis[t] = 0;

//     fronth.update(s, 0);
//     backh.update(t, 0);
//     int ti = 0;
//     int meet = -1;

//     while (fronth.size() || backh.size()) {
//         int element_s, key_s, element_t, key_t;

//         if (fronth.size()) {
//             fronth.extract_min(element_s, key_s);
//             closed_s[element_s] = true;

//             if (backdis[element_s] < 1e9) {
//                 answer = std::min(answer, frontdis[element_s] + backdis[element_s]);
//                 meet = element_s;
//                 break; // Found a meeting point
//             }

//             for (auto edge : graph[element_s]) {
//                 int v = edge.first;
//                 if (CH_order[v] < CH_order[element_s]) continue;
//                 if (closed_t[v]) {
//                     answer = std::min(answer, frontdis[element_s] + edge.second + backdis[v]);
//                     meet = v;
//                     break; // Found a meeting point
//                 }
//                 if (!closed_s[v] && frontdis[v] > frontdis[element_s] + edge.second) {
//                     frontdis[v] = frontdis[element_s] + edge.second;
//                     pre_s[v] = element_s;
//                     stc_s[v] = shortcut_map.find({element_s, v}) != shortcut_map.end();
//                     fronth.update(v, frontdis[v]);
//                 }
//             }
//         }

//         if (backh.size()) {
//             backh.extract_min(element_t, key_t);
//             closed_t[element_t] = true;

//             if (frontdis[element_t] < 1e9) {
//                 answer = std::min(answer, frontdis[element_t] + backdis[element_t]);
//                 meet = element_t;
//                 break; // Found a meeting point
//             }

//             for (auto edge : graph[element_t]) {
//                 int v = edge.first;
//                 if (CH_order[v] < CH_order[element_t]) continue;
//                 if (closed_s[v]) {
//                     answer = std::min(answer, frontdis[v] + edge.second + backdis[element_t]);
//                     meet = v;
//                     break; // Found a meeting point
//                 }
//                 if (!closed_t[v] && backdis[v] > backdis[element_t] + edge.second) {
//                     backdis[v] = backdis[element_t] + edge.second;
//                     pre_t[v] = element_t;
//                     stc_t[v] = shortcut_map.find({element_t, v}) != shortcut_map.end();
//                     backh.update(v, backdis[v]);
//                 }
//             }
//         }
//     }

//     std::vector<std::pair<double, double>> path;
//     std::vector<std::pair<double, double>> search;

//     if (frontdis[t] == 1e9) {
//         return {path, search};
//     }

//     // Trace path from s to meet
//     int current = meet;
//     while (current != s) {
//         if (stc_s[current]) {
//             auto& inside = shortcut_map[{pre_s[current], current}];
//             for (auto it = inside.rbegin(); it != inside.rend(); ++it) {
//                 path.push_back({lat[*it], lon[*it]});
//             }
//         } else {
//             path.push_back({lat[current], lon[current]});
//         }
//         current = pre_s[current];
//     }
//     path.push_back({lat[s], lon[s]});
//     std::reverse(path.begin(), path.end());

//     // Trace path from meet to t
//     current = meet;
//     std::vector<std::pair<double, double>> reverse_path;
//     while (current != t) {
//         if (stc_t[current]) {
//             auto& inside = shortcut_map[{pre_t[current], current}];
//             for (auto it = inside.begin(); it != inside.end(); ++it) {
//                 reverse_path.push_back({lat[*it], lon[*it]});
//             }
//         } else {
//             reverse_path.push_back({lat[current], lon[current]});
//         }
//         current = pre_t[current];
//     }
//     reverse_path.push_back({lat[t], lon[t]});
//     path.insert(path.end(), reverse_path.begin(), reverse_path.end());

//     return {path, search};
// }
