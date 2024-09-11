#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
// #include "CH.cpp"
#include "PLL.cpp"
#include "chrono"
int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    int start, end;
    std::cin >> start >> end;
    
    std::cout << "Start: " << start << ", End: " << end << std::endl;

    // Initialize graph
    Graph myGraph("uploads/nodes.txt", "uploads/graph.txt");
    std::cout<<"ready"<<'\n';
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << " Graph initialized." << std::endl;

    std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> result;
    int length;
   
    // Run dijkstra_way
    

    // Run dijkstra
    // std::cout << "Before Dijkstra length" << std::endl;
    length = myGraph.q2Hop(start, end);
    std::cout << "length completed: " <<length<< '\n';
    
    // std::cout << "Before Dijkstra way" << std::endl;
    result = myGraph.PLL_way(start, end, length);
    
    std::cout << "way completed" << '\n';
    // Extract path and search result
    std::vector<std::pair<double, double>> path = result.first;
    std::vector<std::pair<double, double>> search = result.second;

    // Write path to file
    std::ofstream pathFile("path.txt");
    std::cout << "Writing path to file" << std::endl;
    if (pathFile.is_open()) {
        pathFile << length << "\n";
        for (const auto& p : path) {
            pathFile << p.second << " " << p.first << "\n";
        }
        pathFile.close();
        std::cout << "Path file written." << '\n';
    } else {
        std::cerr << "Error opening path.txt file" << std::endl;
        return 1;
    }

    // Write search to file
    std::ofstream searchFile("search.txt");
    std::cout << "Writing search to file" << std::endl;
    if (searchFile.is_open()) {
        for (const auto& points : search) {
            searchFile << points.second << " " << points.first << "\n";
        }
        searchFile.close();
        std::cout << "Search file written." << '\n';
    } else {
        std::cerr << "Error opening search.txt file" << std::endl;
        return 1;
    }

    std::cout << "Length: " << length << std::endl;

    return 0;
}






