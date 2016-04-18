//
// Created by Seth Wiesman
//

#ifndef METRICSPACEDATASTRUCTURES_METRICTREE_HPP
#define METRICSPACEDATASTRUCTURES_METRICTREE_HPP

#include <vector>
#include <algorithm>
#include <memory>

template<typename T, double(*distance)(const T &, const T &)>
class MetricTree  {
    struct Node {
        T point;

        double innerRadius;
        double outerRadius;
        double boundingRadius;

        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;

        Node(T point) :
                point(point),
                innerRadius(0),
                outerRadius(0),
                boundingRadius(0) { }

        void collect(std::vector<T> &result);
        void search(
                const T &target,
                const double radius,
                std::vector<T> &result
        );
    };

    using node_itr = typename std::vector<std::shared_ptr<typename MetricTree<T, distance>::Node>>::iterator;

    std::shared_ptr<Node> root;

    std::shared_ptr<Node> buildTree(const node_itr begin, const node_itr end);

public:
    MetricTree(std::vector<T> points);
    std::vector<T> search(const T &target, const double radius);
};

template<typename T, double(*distance)(const T &, const T &)>
MetricTree<T, distance>::MetricTree(std::vector<T> points) {
    std::vector<std::shared_ptr<Node>> nodes;
    nodes.reserve(points.size());

    for (auto &point : points) {
        nodes.push_back(std::make_shared<Node>(point));
    }
    this->root = buildTree(nodes.begin(), nodes.end());
}

template<typename T, double(*distance)(const T &, const T &)>
std::vector<T> MetricTree<T, distance>::search(const T &target, const double radius) {

    std::vector<T> result;
    if (this->root != nullptr) {
        this->root->search(target, radius, result);
    }
    return result;
}

template<typename T, double(*distance)(const T &, const T &)>
std::shared_ptr<typename MetricTree<T, distance>::Node>
MetricTree<T, distance>::buildTree(const node_itr begin, const node_itr end) {
    if (begin == end) {
        return nullptr;
    }

    if ((end - begin) == 1) {
        return *begin;
    }

    for (auto itr = begin + 1; itr != end; itr++) {
        (*itr)->innerRadius = distance((*begin)->point, (*itr)->point);
    }

    const auto median = begin + (end - begin) / 2;

    const auto cmp = [](const auto left, const auto right) {
        return left->innerRadius < right->innerRadius;
    };

    std::nth_element(begin + 1, median, end, cmp);

    (*begin)->outerRadius = (*median)->innerRadius;

    const auto pointOnInnerRadius = std::max_element(begin + 1, median, cmp);
    const auto pointOnOuterBound  = std::max_element(median, end, cmp);

    (*begin)->innerRadius    = (*pointOnInnerRadius)->innerRadius;
    (*begin)->boundingRadius = (*pointOnOuterBound)->innerRadius;

    (*begin)->left  = buildTree(begin + 1, median);
    (*begin)->right = buildTree(median, end);

    return *begin;
}

template<typename T, double(*distance)(const T &, const T &)>
void MetricTree<T, distance>::Node::search(const T &target, const double radius, std::vector<T> &result) {

    const auto dist = distance(this->point, target);

    if (dist <= radius) {
        result.push_back(this->point);
    }

    if (this->left != nullptr) {
        if (dist + this->innerRadius <= radius) {
            this->left->collect(result);
        } else if (dist - radius <= this->innerRadius) {
            this->left->search(target, radius, result);
        }
    }

    if (this->right != nullptr) {
        if (dist + this->boundingRadius <= radius) {
            this->right->collect(result);
        } else if (dist + radius >= this->outerRadius) {
            this->right->search(target, radius, result);
        }
    }
}

template<typename T, double(*distance)(const T &, const T &)>
void MetricTree<T, distance>::Node::collect(std::vector<T>& result) {
    result.push_back(this->point);
    if (this->left != nullptr) {
        this->left->collect(result);
    }
    if (this->right != nullptr) {
        this->right->collect(result);
    }
}

#endif //METRICSPACEDATASTRUCTURES_METRICTREE_HPP
