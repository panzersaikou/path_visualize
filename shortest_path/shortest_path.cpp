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
    std::cout<<"!";
    //ALT_preprocess
    Preprocessed_path.resize(20);
    for(int i=0;i<20;i++)
    {
        preprocess(i,number_of_nodes/20*i);
    }
    std::cout<<real_dis(lat[0],lon[0],lat[1],lon[1])<<std::endl;
    clock_t start=clock();

    //preprocess2Hop();
    H2H_preprocess();
    
    //CH_preprocess();
    std::cout<<ALT(0,2)<<std::endl;
    std::cout<<Astar(0,2)<<std::endl;
    std::cout<<dijkstra(0,2)<<std::endl;
    std::cout<<Bidijkstra(0,2)<<std::endl;
    std::cout<<H2H(0,2)<<std::endl;
    //self_test();
    //check();
    //makequery();
    work();
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
        for (const auto& neighbora : graph[node])
        {
            int v=neighbora.first;
            if(contracted_nodes[v])continue;
            for (const auto& neighborb : graph[node])
            {
                int w=neighborb.first;
                int shortcut=neighbora.second+neighborb.second;
                if(contracted_nodes[w])continue;
                if(v==w)continue;
                
                if(exist_line[v][w])
                {
                    for(auto& edge:graph[v])
                    {
                        if(edge.first==w)
                        {
                            edge.second=std::min(edge.second,shortcut);
                            break;
                        }
                    }
                    for(auto& edge:graph[w])
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
                    graph[v].push_back(std::make_pair(w,shortcut));
                    graph[w].push_back(std::make_pair(v,shortcut));
                    nodes[v].degree++;
                    nodes[w].degree++;
                    exist_line[v][w]=exist_line[w][v]=1;
                }
            }
        }
        for (const auto& neighbora : graph[node])
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
        for(auto& edge:graph[i])
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

void Graph::preprocess2Hop()
{
    clock_t start=clock();
    sort(nodes.begin(),nodes.end(),comparenodes);
    for(int i=0;i < number_of_nodes; i++)
    {
        order[nodes[i].osmid]=i;
    }
    for (int i = 0; i < number_of_nodes; ++i) {
        prunedBFS(nodes[i].osmid);
        if(i%(number_of_nodes/1000)==0)
        {
            std::cout<<nodes[i].degree<<" "<<nodes[i].osmid<<std::endl;
            clock_t end = clock();
            double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
            std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        }
    }
}

void Graph::CH_preprocess()
{
    heap q(number_of_nodes);
    for (int i = 0; i < number_of_nodes; i++) {
        q.update(i, nodes[i].degree);
    }
    clock_t start=clock();
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
        for (const auto& neighbora : graph[node])
        {
            int v=neighbora.first;
            if(contracted_nodes[v])continue;
            for (const auto& neighborb : graph[node])
            {
                int w=neighborb.first;
                int shortcut=neighbora.second+neighborb.second;
                if(contracted_nodes[w])continue;
                if(v==w)continue;
                if(exist_line[v][w])
                {
                    for(auto& edge:graph[v])
                    {
                        if(edge.first==w)
                        {
                            edge.second=std::min(edge.second,shortcut);
                            break;
                        }
                    }
                    for(auto& edge:graph[w])
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
                    graph[v].push_back(std::make_pair(w,shortcut));
                    graph[w].push_back(std::make_pair(v,shortcut));
                    nodes[v].degree++;
                    nodes[w].degree++;
                    exist_line[v][w]=exist_line[w][v]=1;
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

        if (cnt % (number_of_nodes / 100) == 0) {
            std::cout << priority << " " << node <<" "<<cnt<< std::endl;
            clock_t end = clock();
            double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
            std::cout << "Elapsed CPU time: " << elapsed_secs << " seconds." << std::endl;
        }
    }
    std::cout<<"done preprocess"<<std::endl;
    for(int i=0;i<number_of_nodes;i++)contracted_nodes[i]=0;
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
void Graph::check()
{
    srand(time(0));
    for(int i=0;i<1000;i++)
    {
        int s=rand()%number_of_nodes,t=rand()%number_of_nodes;
        while(s==t)t=rand()%number_of_nodes;
        //s=11115,t=32409;
        int d1=ALT(s,t),d2=Astar(s,t),d3=dijkstra(s,t),d4=Bidijkstra(s,t);
        //std::cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<std::endl;
        /*std::cout<<s<<" "<<t<<std::endl;
        int x=t;
        while(x!=s)
        {
            std::cout<<x<<" -> ";
            x=pre[x];
        }
        std::cout<<s<<std::endl;*/
        assert(d1==d2);
        assert(d2==d3);
        assert(d3==d4);
    }
}

void Graph::makequery()
{
    double min_lat,min_lon,max_lat,max_lon;
    min_lat=max_lat=lat[0];
    min_lon=max_lon=lon[0];
    for(int i=0;i<number_of_nodes;i++)
    {
        min_lat=std::min(min_lat,lat[i]);
        min_lon=std::min(min_lon,lon[i]);
        max_lat=std::max(max_lat,lat[i]);
        max_lon=std::max(max_lon,lon[i]);
    }
    double lmin=1000,lmax=real_dis(min_lat,min_lon,max_lat,max_lon);
    std::cout<<lmin<<" "<<lmax<<std::endl;
    double x=std::pow(lmax/lmin,0.1);
    int tot=0;
    queries.resize(10);
    while(tot<1000)
    {
        if(tot%100==0)
        {
            for(int i=0;i<10;i++)std::cout<<queries[i].size()<<" ";
            std::cout<<std::endl;
        }
        int s=rand()*rand()%number_of_nodes,t=rand()*rand()%number_of_nodes;
        while(s==t)t=rand()*rand()%number_of_nodes;
        int dis=real_dis(lat[s],lon[s],lat[t],lon[t]);
        for(int i=1;i<=10;i++)
        {
            int left=(int)(lmin*pow(x,i-1)),right=(int)(lmin*pow(x,i));
            if(dis>left&&dis<=right)
            {
                if(queries[i-1].size()==100)break;
                queries[i-1].push_back(std::make_pair(s,t));
                tot++;
                break;
            }
        }
    }
    std::ofstream file("GZ.query");
    file<<1000<<std::endl;
    for(int i=0;i<10;i++)
    {
        for(int j=0;j<100;j++)
        {
            file<<queries[i][j].first<<" "<<queries[i][j].second<<std::endl;
            if(j==0)std::cout<<Astar(queries[i][j].first,queries[i][j].second)<<std::endl;
        }
    }
    file.close();
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
    std::vector<int>vis;
    vis.resize(number_of_nodes,0);
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
                pre[v]=element;
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
        //cnt++;
        int element,key;
        q.extract_min(element,key);
        if(element==t)break;
        vis[element]=1;
        //std::cout<<element<<std::endl;
        for(auto edge : graph[element])
        {
            int v=edge.first;
            if(vis[v])continue;
            if(distance[v]>distance[element]+edge.second)
            {
                distance[v]=distance[element]+edge.second;
                pre[v]=element;
                //std::cout<<v<<std::endl;
                //q.update(v,distance[v]+query_h(v,t));
                q.update(v,distance[v]);
                //std::cout<<v<<std::endl;
                vis[v]=0;
            }
        }
    }
    //std::cout<<t<<std::endl;
    if(t==number_of_nodes)return 0;
    return distance[t];
}

int Graph::Bidijkstra(int s,int t)
{
    //std::cout<<s<<" "<<t<<std::endl;
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
    //using dijistra and A*
    fronth.update(s,0);
    backh.update(t,0);
    while(fronth.size()&&backh.size())
    {
        //if(fronth.top()+backh.top()>=answer)break;
        cnt+=2;
        int element_s,key_s,element_t,key_t;
        fronth.extract_min(element_s,key_s);
        backh.extract_min(element_t,key_t);
        closed_s[element_s]=true;
        if(closed_t[element_s])break;
        //answer=std::min(answer,frontdis[element_s]+backdis[element_s]);
        //std::cout<<element<<std::endl;
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
        //answer=std::min(answer,frontdis[element_t]+backdis[element_t]);
        //std::cout<<element<<std::endl;

        closed_t[element_t]=true;
        if(closed_s[element_t])break;
        for(auto edge : graph[element_t])
        {
            int v=edge.first;
            //std::cout<<v<<std::endl;
            if(closed_s[v])answer=std::min(answer,frontdis[v]+edge.second+backdis[element_t]);
            if(!closed_t[v]&& backdis[v]>backdis[element_t]+edge.second)
            {
                backdis[v]=backdis[element_t]+edge.second;
                backh.update(v,backdis[v]);
                //std::cout<<v<<std::endl;
            }
        }
    }
    //std::cout<<t<<std::endl;
    return answer;

}

int Graph::CH(int s,int t)
{
    if(s==t)return 0;
    //std::cout<<s<<" "<<t<<std::endl;
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
    //using dijistra and A*
    fronth.update(s,0);
    backh.update(t,0);
    int ti=0;
    while(fronth.size()||backh.size())
    {
        //if(fronth.top()+backh.top()>=answer)break;
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
                if(closed_t[v])answer=std::min(answer,frontdis[element_s]+edge.second+backdis[v]);
                if(!closed_s[v]&& frontdis[v]>frontdis[element_s]+edge.second)
                {
                    frontdis[v]=frontdis[element_s]+edge.second;
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

const double EARTH_RADIUS_KM = 6371.0;
int Graph::real_dis(double lat1, double lon1, double lat2, double lon2)
{
    //std::cout<<lat1<<" "<<lon1<<" "<<lat2<<" "<<lon2<<std::endl;
    //return 0;
    /*int lat=(int)(abs(lat1-lat2)*111319);
	int lon=(int)(abs(lon1-lon2)*102175);
	int min,max;
	min=(lat>lon)?lon:lat;
	max=(lat>lon)?lat:lon;
	int approx=max*1007+min*441;
	if(max<(min<<4))
		approx-=max*40;
	return ((approx+512)>>10);*/


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
                pre[v]=x;
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