
// #include "0-index.cpp"

#include "CH.cpp"
int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <start_node> <end_node>" << std::endl;
        return 1;
    }
    int start = std::stoi(argv[1]);
    int end = std::stoi(argv[2]);
    std::string Filepath_G = argv[3];
    std::string Filepath_N = argv[4];
    clock_t t1 = clock();
    Graph myGraph(Filepath_N,Filepath_G);
    clock_t t2 = clock();
    std::pair<std::vector<std::pair<double,double>>,std::vector<std::pair<double, double>>> result;
    int length;
    clock_t t5 = clock();
    length = myGraph.CH(start,end);
    clock_t t6 = clock();
    clock_t t3 = clock();
    result = myGraph.CH_way(start,end);
    clock_t t4 = clock();

    
    std::vector<std::pair<double,double>>path = result.first;
    std::vector<std::pair<double,double>>search = result.second;
    std::ofstream pathFile("path.txt");
    

    if (pathFile.is_open()) {
        pathFile<<length<<"\n";
        for (const auto& p : path) {
            pathFile << p.second << " " << p.first << "\n";
        }
        pathFile.close();
    } else {
        std::cerr << "Error opening file" << std::endl;
        
        return 1;
    }

    std::ofstream searchFile("search.txt");
    if (searchFile.is_open()) {
        for (const auto& points : search) {
            searchFile << points.second << " " << points.first << "\n";
        }
        searchFile.close();
    } else {
        std::cerr << "Error opening search.txt file" << std::endl;
        
        return 1;
    }

    std::cout<<"time for prepare:"<<static_cast<double> (t2-t1)/CLOCKS_PER_SEC<<"seconds"<<'\n'<<"time for find distance:"<<static_cast<double>(t6-t5)/CLOCKS_PER_SEC<<"seconds"<<'\n'<<"time for construct way:"<<static_cast<double>(t4-t3)/CLOCKS_PER_SEC<<"microsecond.";

    return 0;
}

