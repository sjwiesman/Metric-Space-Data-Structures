#include <iostream>
#include <random>
#include "MetricTree.hpp"
#include "FaSST.hpp"
#include "FaSSTGating.hpp"

using Point = std::vector<int>;

std::vector<Point> generatePoints(std::vector<int>::size_type size, std::vector<int>::size_type dimensionality) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,100);

    std::vector<Point> points(size, Point(dimensionality, 0));
    for (auto i = 0; i < size; i++) {
        for (auto j = 0; j < dimensionality; j++) {
            points[i][j] = distribution(generator);
        }
    }

    return points;
}

double euclideanDistance(const Point& p1, const Point& p2) {
    auto sum = 0.0;
    for (auto i = 0; i < p1.size(); i++) {
        sum += (p1[i] + p2[i])*(p1[i] + p2[i]);
    }

    return sum / p1.size();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Improper Usage" << std::endl;
        std::cout << argv[0] << " <num points> <dimensionality>" << std::endl;
    }

    Point::size_type size, dimensionality;
    sscanf(argv[1], "%zu", &size);
    sscanf(argv[2], "%zu", &dimensionality);

    auto radius = 5;
    auto origin = Point(dimensionality, 0);
    auto points = generatePoints(size, dimensionality);

    MetricTree<Point, euclideanDistance>  metricTree(points);
    auto metricTreeResultSet  = metricTree.search(origin, radius);

    while (metricTreeResultSet.size() == 0) {
        radius *= 1.5;
        metricTreeResultSet = metricTree.search(origin, radius);
    }

    std::cout << "Radius" << radius << std::endl << std::endl;

    std::cout << "Metric Tree Result" << std::endl;
    for (auto& point : metricTreeResultSet) {
        std::cout << '(';
        for(auto& elem : point) {
            std:: cout << elem << ',';
        }
        std::cout << "\b \b" << ')' << std::endl;
    }

    FaSST<Point, euclideanDistance> faSST(points);
    auto faSSTResultSet  = faSST.search(origin, radius);
    std::cout << "Fast Similarity Search Tree Result" << std::endl;
    for (auto& point : faSSTResultSet) {
        std::cout << '(';
        for(auto& elem : point) {
            std:: cout << elem << ',';
        }
        std::cout << "\b \b" << ')' << std::endl;
    }

    FaSSTGating<Point, euclideanDistance> faSSTGating(points);
    auto faSSTGatingResultSet = faSSTGating.search(origin, radius);
    std::cout << "Fast Similarity Search Tree with Gating Result" << std::endl;
    for (auto& point : faSSTGatingResultSet) {
        std::cout << '(';
        for(auto& elem : point) {
            std:: cout << elem << ',';
        }
        std::cout << "\b \b" << ')' << std::endl;
    }

    return 0;
}