//
// Created by Seth Wiesman on 4/18/16.
//

#ifndef METRICSPACEDATASTRUCTURES_FASSTGATING_HPP
#define METRICSPACEDATASTRUCTURES_FASSTGATING_HPP

#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>

template<typename T, double (*distance)(const T &, const T &)>
class FaSSTGating {

    struct Annulus {
        double shortRadius;
        double longRadius;

        Annulus(double shortRadius, double longRadius) :
                shortRadius(shortRadius),
                longRadius(longRadius) { }
    };

    struct Node {
        T point;

        std::vector<Annulus> leftAnnuli;
        std::vector<Annulus> rightAnnuli;

        std::vector<double> pivots;

        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;

        Node(T point) : point(point) { }

        void collect(std::vector<T>& result);
        void search(
                const T &target,
                const double radius,
                std::vector<double> &pivots,
                std::vector<T> &result
        );

        bool collectLeft(std::vector<double>& pivots, double radius);
        bool goLeft(std::vector<double> &pivots, double radius);

        bool collectRight(std::vector<double>& pivots, double radius);
        bool goRight(std::vector<double> &pivots, double radius);
    };

    using node_itr = typename std::vector<std::shared_ptr<typename FaSSTGating<T, distance>::Node>>::iterator;

    std::shared_ptr<Node> root;

    std::shared_ptr<Node> buildTree(const node_itr begin, const node_itr high);

public:
    FaSSTGating(std::vector<T> points);
    std::vector<T> search(const T &target, const double radius);
};

template<typename T, double (*distance)(const T &, const T &)>
FaSSTGating<T, distance>::FaSSTGating(std::vector<T> points) {
    std::vector<std::shared_ptr<Node>> nodes;
    nodes.reserve(points.size());

    for (auto &point : points) {
        nodes.push_back(std::make_shared<Node>(point));
    }

    this->root = buildTree(nodes.begin(), nodes.end());
}

template<typename T, double (*distance)(const T &, const T &)>
std::vector<T> FaSSTGating<T, distance>::search(const T &target, double radius) {
    std::vector<T> result;
    std::vector<double> pivots;

    if (this->root != nullptr) {
        this->root->search(target, radius, pivots, result);
    }

    return result;
}

template<typename T, double (*distance)(const T &, const T &)>
std::shared_ptr<typename FaSSTGating<T, distance>::Node>
FaSSTGating<T, distance>::buildTree(const node_itr begin, const node_itr end) {
    if (begin == end) {
        return nullptr;
    }

    if ((end - begin) == 1) {
        return *begin;
    }

    for (auto itr = begin + 1; itr != end; itr++) {
        auto dist = distance((*begin)->point, (*itr)->point);
        (*itr)->pivots.push_back(dist);
    }

    const auto median = begin + (end - begin) / 2;

    std::nth_element(begin + 1, median, end, [](auto left, auto right) {
        return left->pivots.back() < right->pivots.back();
    });

    (*begin)->left  = buildTree(begin + 1, median);
    (*begin)->right = buildTree(median, end);

    unsigned long depth = 0;
    if ((*begin)->left != nullptr) {
        depth = (*begin)->left->pivots.size();
    } else if ((*begin)->right != nullptr) {
        depth = (*begin)->right->pivots.size();
    }

    for (auto i = 0; i < depth; i++) {
        const auto cmp = [i](auto left, auto right) {
            return left->pivots[i] < right->pivots[i];
        };

        if ((*begin)->left != nullptr) {
            const auto nearest  = std::min_element(begin + 1, median, cmp);
            const auto furthest = std::max_element(begin + 1, median, cmp);

            const Annulus annulus((*nearest)->pivots[i], (*furthest)->pivots[i]);
            (*begin)->leftAnnuli.push_back(annulus);
        }

        if ((*begin)->right != nullptr) {
            const auto nearest  = std::min_element(median, end, cmp);
            const auto furthest = std::max_element(median, end, cmp);

            const Annulus annulus((*nearest)->pivots[i], (*furthest)->pivots[i]);
            (*begin)->rightAnnuli.push_back(annulus);
        }
    }
    return *begin;
}

template<typename T, double(*distance)(const T &, const T &)>
void FaSSTGating<T, distance>::Node::search(
        const T &target,
        const double radius,
        std::vector<double> &pivots,
        std::vector<T> &result
) {
    auto min = 0.0;
    auto max = std::numeric_limits<double>::max();

    for (auto i = 0; i < this->pivots.size(); i++) {
        if (pivots[i] == std::numeric_limits<double>::max()) {
            continue;
        }

        const auto minRange = std::fabs(this->pivots[i] - pivots[i]);
        if (min < minRange) {
            min = minRange;
        }

        const auto maxRange = this->pivots[i] + pivots[i];
        if (maxRange < max) {
            max = maxRange;
        }
    }

    if (max <= radius) {
        result.push_back(this->point);
        pivots.push_back(std::numeric_limits<double>::max());
    } else if (radius < min) {
        pivots.push_back(std::numeric_limits<double>::max());
    } else {
        const auto dist = distance(target, this->point);
        if (dist <= radius) {
            result.push_back(this->point);
        }

        pivots.push_back(dist);
    }

    if (this->left != nullptr) {
        if (this->collectLeft(pivots, radius)) {
            this->left->collect(result);
        } else if (this->goLeft(pivots, radius)) {
            this->left->search(target, radius, pivots, result);
        }
    }

    if (this->right != nullptr) {
        if (this->collectRight(pivots, radius)) {
            this->right->collect(result);
        } else if (this->goRight(pivots, radius)) {
            this->right->search(target, radius, pivots, result);
        }
    }
}

template<typename T, double(*distance)(const T&, const T&)>
bool FaSSTGating<T, distance>::Node::collectLeft(
        std::vector<double>& pivots,
        double radius
) {
    for (auto i = 0; i < this->leftAnnuli.size(); i++) {
        if (pivots[i] == std::numeric_limits<double>::max()) {
            continue;
        }
        if (pivots[i] + this->leftAnnuli[i].longRadius <= radius) {
            return true;
        }
    }
    return false;
}

template<typename T, double(*distance)(const T&, const T&)>
bool FaSSTGating<T, distance>::Node::goLeft(std::vector<double> &pivots, double radius) {
    for (auto i = 0; i < this->pivots.size(); i++) {
        if (pivots[i] == std::numeric_limits<double>::max()) {
            continue;
        }
        if (pivots[i] < this->leftAnnuli[i].shortRadius) {
            if (pivots[i] + radius < this->leftAnnuli[i].shortRadius) {
                return false;
            }
        } else if (pivots[i] - radius > this->leftAnnuli[i].longRadius) {
            if (this->leftAnnuli[i].longRadius < pivots[i]) {
                return false;
            }
        }
    }

    return true;
}

template<typename T, double(*distance)(const T&, const T&)>
bool FaSSTGating<T, distance>::Node::collectRight(std::vector<double>& pivots, double radius) {
    for (auto i = 0; i < this->rightAnnuli.size(); i++) {
        if (pivots[i] == std::numeric_limits<double>::max()) {
            continue;
        }

        if (pivots[i] + this->rightAnnuli[i].longRadius <= radius) {
            return true;
        }
    }
    return false;
}

template<typename T, double(*distance)(const T&, const T&)>
bool FaSSTGating<T, distance>::Node::goRight(
        std::vector<double> &pivots,
        double radius
) {
    for (auto i = 0; i < this->pivots.size(); i++) {
        if (pivots[i] == std::numeric_limits<double>::max()) {
            continue;
        }
        if (pivots[i] < this->rightAnnuli[i].shortRadius) {
            if (pivots[i] + radius < this->rightAnnuli[i].shortRadius) {
                return false;
            }
        } else if (pivots[i] - radius > this->rightAnnuli[i].longRadius) {
            if (this->rightAnnuli[i].longRadius < pivots[i]) {
                return false;
            }
        }
    }

    return true;
}

template<typename T, double(*distance)(const T&, const T&)>
void FaSSTGating<T, distance>::Node::collect(std::vector<T>& result) {
    result.push_back(this->point);
    if (this->left != nullptr) {
        this->left->collect(result);
    }

    if (this->right != nullptr) {
        this->right->collect(result);
    }
}

#endif //METRICSPACEDATASTRUCTURES_FASSTGATING_HPP
