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

Graph::Graph() 
{
    readin();
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

void Graph::readin() 
{
    int bo=1;
    exist_line.clear();
    //readin a traffic graph
    std::ifstream file("GZ.co");
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
        H2H_graph.resize(number_of_nodes);
        cnt_of_path.resize(number_of_nodes);
        for (int i = 0; i < number_of_nodes; i++) {
            file >> ix >> dx >> dy;
            lon[i] = dx/1000000.0;
            lat[i] = dy/1000000.0;
            nodes[i].osmid=i;
        }
        file.close();
    } else {
        bo=0;
        std::cout << "Unable to open file." << std::endl;
    }
    if(!bo)return;

    //get the information of streets
    file.open("GZ");

    if (file.is_open()) {
        file >> number_of_nodes>> number_of_line;
        for (int i = 0; i < number_of_line; i++) {
            file>>ix>>iy>>dx;
            if(exist_line[ix][iy])
            {
                for(auto& edge:H2H_graph[ix])
                {
                    if(edge.first==iy)
                    {
                        edge.second=std::min(edge.second,(int)dx);
                        break;
                    }
                }
                for(auto& edge:H2H_graph[iy])
                {
                    if(edge.first==ix)
                    {
                        edge.second=std::min(edge.second,(int)dx);
                        break;
                    }
                }
            }
            else
            {
                H2H_graph[ix].push_back(std::make_pair(iy,dx));
                H2H_graph[iy].push_back(std::make_pair(ix,dx));
            }
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
    for(int i=0;i<number_of_nodes;i++)
    {
        if(graph[i]!=H2H_graph[i])
        {
            for(auto poi:graph[i])std::cout<<poi.first<<" "<<poi.second<<" ";
            std::cout<<std::endl;
            for(auto poi:H2H_graph[i])std::cout<<poi.first<<" "<<poi.second<<" ";
            std::cout<<std::endl;
            return;
        }
    }
    std::cout<<"!";
    H2H_preprocess();
    
    work();
    self_test();
}

void Graph::work()
{
    cnt=0;
    std::ifstream infile("GZ.query");
    infile>>number_of_query;
    test.resize(number_of_query);
    std::cout<<number_of_query<<std::endl;
    std::ofstream outfile("GZ.answerH2H");
    clock_t start = clock();
    for(int i=0;i<number_of_query;i++)
    {
        int s,t;
        infile>>s>>t;
        outfile<<H2H(s,t)<<std::endl;
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

void Graph::H2H_preprocess()
{
    lgsiz=log(number_of_nodes)/log(2)+2;
    for(int i=0;i<number_of_nodes;i++)fa[i].resize(lgsiz);
    std::stack<int>stk;
    //make a DPtree
    heap q(number_of_nodes);
    for (int i = 0; i < number_of_nodes; i++) {
        q.update(i, nodes[i].degree);
    }
    clock_t start=clock();
    while (!q.empty()) {
        int node, priority;
        q.extract_min(node, priority);
        stk.push(node);
        H2H_order[node]=cnt;
        cnt++;
        contracted_nodes[node]=1;
        //std::cout<<"test"<<std::endl;
        for (const auto& neighbora : H2H_graph[node])
        {
            int v=neighbora.first;
            if(contracted_nodes[v])continue;
            for (const auto& neighborb : H2H_graph[node])
            {
                int w=neighborb.first;
                int shortcut=neighbora.second+neighborb.second;
                if(contracted_nodes[w])continue;
                if(v==w)continue;
                
                if(exist_line[v][w])
                {
                    for(auto& edge:H2H_graph[v])
                    {
                        if(edge.first==w)
                        {
                            edge.second=std::min(edge.second,shortcut);
                            break;
                        }
                    }
                    for(auto& edge:H2H_graph[w])
                    {
                        if(edge.first==v)
                        {
                            edge.second=std::min(edge.second,shortcut);
                            break;
                        }
                    }
                }
                else
                {
                    H2H_graph[v].push_back(std::make_pair(w,shortcut));
                    H2H_graph[w].push_back(std::make_pair(v,shortcut));
                    nodes[v].degree++;
                    nodes[w].degree++;
                    exist_line[v][w]=exist_line[w][v]=1;
                }
            }
        }
        for (const auto& neighbora : H2H_graph[node])
        {
            int v=neighbora.first;
            if(contracted_nodes[v])continue;
            X[node].push_back(v);
            nodes[v].degree--;
            q.update(v,nodes[v].degree);
        }
        X[node].push_back(node);

        if (cnt % (number_of_nodes / 100) == 0) {
            std::cout << priority << " " << node <<" "<<cnt<< std::endl;
            clock_t end = clock();
            double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
            std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        }
    }
    for(int i=0;i<number_of_nodes;i++)
    {
        //sort(X[i].begin(),X[i].end(),compareX);
        sort(X[i].begin(),X[i].end(),[this](int a, int b) { return this->compareX(a, b); });
        //std::cout<<"test"<<std::endl;
        fi_dis[i].resize(X[i].size());
        if(X[i].size()==1)
        {
            fa[i][0]=-1;
            continue;
        }
        fa[i][0]=X[i][X[i].size()-2];
        for(auto& edge:H2H_graph[i])
        {
            distance[edge.first]=edge.second;
        }
        for(int j=0;j<X[i].size()-1;j++)
        {
            fi_dis[i][j]=distance[X[i][j]];
        }
        fi_dis[i][X[i].size()-1]=0;
    }
    //get pos&dis
    std::vector<int>anc_pos;
    anc_pos.resize(number_of_nodes);
    int tot=0;
    while(stk.size())
    {
        int node=stk.top();stk.pop();
        if(fa[node][0]!=-1)
        {
            //std::cout<<"hasfa"<<std::endl;
            dep[node]=dep[fa[node][0]]+1;
            for(int i=0;i<ancestor[fa[node][0]].size();i++)
            {
                int v=ancestor[fa[node][0]][i];
                ancestor[node].push_back(v);
                anc_pos[v]=i;
            }
        }
        ancestor[node].push_back(node);
        anc_pos[node]=ancestor[node].size()-1;
        pos[node].resize(X[node].size());
        //dis[node].resize(ancestor[node].size());
        for(int i=0;i<X[node].size();i++)
        {
            pos[node][i]=anc_pos[X[node][i]];
        }
        //std::cout<<H2H_order[node]<<std::endl;
        /*for(int i=0;i<X[node].size()-1;i++)
        {
            std::cout<<pos[node][i]<<std::endl;
        }
        std::cout<<ancestor[node].size()<<std::endl;
        std::cout<<"test"<<std::endl;*/
        for(int i=0;i<ancestor[node].size()-1;i++)
        {
            int mn_dis=1e9;
            for(int j=0;j<X[node].size()-1;j++)
            {
                int d;
                if(pos[node][j]>i)d=dis[X[node][j]][i];
                else d=dis[ancestor[node][i]][pos[node][j]];
                //std::cout<<"?"<<std::endl;
                mn_dis=std::min(mn_dis,d+fi_dis[node][j]);
                //std::cout<<"!"<<std::endl;
            }
            dis[node].push_back(mn_dis);
        }
        //std::cout<<++tot<<std::endl;
        dis[node].push_back(0);
    }
    //std::cout<<"test"<<std::endl;
    for(int i=1;i<lgsiz;i++)
    {
        for(int j=0;j<number_of_nodes;j++)
        {
            if(fa[j][i-1]==-1)fa[j][i]=-1;
            else fa[j][i]=fa[fa[j][i-1]][i-1];
        }
    }
    //std::cout<<"test"<<std::endl;
}

int Graph::LCA(int s,int t)
{
    if(dep[s]<dep[t])std::swap(s,t);
    for(int j=lgsiz-1;j>=0;j--)
    {
        if(fa[s][j]!=-1&&dep[fa[s][j]]>=dep[t])s=fa[s][j];
    }
    if(s==t)return s;
    for(int j=lgsiz-1;j>=0;j--)
    {
        if(fa[s][j]!=fa[t][j])s=fa[s][j],t=fa[t][j];
    }
    return fa[s][0];
}

int Graph::H2H(int s,int t)
{
    //std::cout<<s<<" "<<t<<std::endl;
    int z=LCA(s,t);
    int ans=1e9;
    for(int i:pos[z])
    {
        ans=std::min(ans,dis[s][i]+dis[t][i]);
    }
    return ans;
}

std::vector<int> Graph::H2H_path(int s,int t,int dist)
{
    int i=0;
    std::vector<int>le,ri,path;
    for(;i<=ancestor[s].size();i++)
    {
        if(dis[s][i]+dis[t][i]==dist)
        {
            break;
        }
    }
    int u=ancestor[s][i];
    //std::cout<<u<<std::endl;
    while(s!=u)
    {
        le.push_back(s);
        for(auto edge:graph[s])
        {
            int v=edge.first;
            if(dep[u]<=dep[v]&&edge.second+dis[v][i]==dis[s][i])
            {
                s=v;
                break;
            }
        }
    }
    //std::cout<<"S"<<std::endl;
    while(t!=u)
    {
        ri.push_back(t);
        for(auto edge:graph[t])
        {
            int v=edge.first;
            if(dep[u]<=dep[v]&&edge.second+dis[v][i]==dis[t][i])
            {
                t=v;
                break;
            }
        }
    }
    //std::cout<<"T"<<std::endl;
    for(auto poin:le)path.push_back(poin);
    path.push_back(u);
    for(i=ri.size()-1;i>=0;i--)path.push_back(ri[i]);
    return path;
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
    //cnt_of_path[s]=1;
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
                //cnt_of_path[v]=cnt_of_path[element];
            }
            /*else if(distance[v]==distance[element]+edge.second)
            {
                cnt_of_path[v]+=cnt_of_path[element];
            }*/
        }
    }
    if(t==number_of_nodes)return 0;
    return distance[t];
}

std::vector<int> Graph::dijkstra_path(int s,int t)
{
    std::stack<int>sta;
    sta.push(t);
    while(t!=s)
    {
        t=pre[t];
        sta.push(t);
    }
    std::vector<int>ve;
    while(sta.size())
    {
        int x=sta.top();sta.pop();
        ve.push_back(x);
    }
    return ve;
}

void Graph::self_test()
{
    for(int i=0;i<1000;i++)
    {
        int s=rand()*rand()%number_of_nodes,t=rand()*rand()%number_of_nodes;
        while(s==t)t=rand()*rand()%number_of_nodes;
        int dis;
        std::vector<int>dij,H2;
        dis=dijkstra(s,t);
        dij=dijkstra_path(s,t);
        dis=H2H(s,t);
        //std::cout<<dis<<"H2Hdis"<<std::endl;
        H2=H2H_path(s,t,dis);
        
        /*if(dij!=H2)
        {
            std::cout<<"wrong"<<std::endl;
            for(auto x:H2)std::cout<<x<<" ";
            std::cout<<std::endl;
            for(auto x:dij)std::cout<<x<<" ";
            std::cout<<std::endl;
            //return;
        }
        else std::cout<<"right"<<std::endl;*/
        
    }
}
