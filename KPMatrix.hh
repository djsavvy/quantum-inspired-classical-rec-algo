#pragma once

#include <cassert>
#include <cmath>
#include <exception>
#include <random>
#include <stdexcept>


class KPVector {
    public:
        KPVector(int n);
        ~KPVector();

        double get(int i);
        void set(int i, double v);
        double squared_norm();
        double sample();

    private:
        int dim;

        struct Node {
            double weight;
            int sign; // 0 if not a leaf, -1 or 1 otherwise

            Node* parent;
            Node* left;
            Node* right;

            Node(Node* p);
            Node(Node* p, double w);
            ~Node();
        };

        Node* find(int i);

        Node* root;
};


KPVector::KPVector(int n) {
    assert(n > 0);

    dim = n;
    root = new KPVector::Node(nullptr);
}

KPVector::Node::Node(Node* p, double w) {
    parent = p;
    weight = w;
    sign = 0;
    left = nullptr;
    right = nullptr;
}

KPVector::Node::Node(Node* p) 
    : KPVector::Node::Node(p, 0) {}

KPVector::Node::~Node() {
    delete left;
    delete right;
}

KPVector::~KPVector() {
    delete root;
}


// Finds pointer to the leaf node for a given index
// If it doesn't yet exist, creates it
KPVector::Node* KPVector::find(int i) {
    // bounds check
    assert(0 <= i && i < dim);

    Node* cur_node = root;
    int low = 0, high = dim;

    // While range has more than one leaf in it, descend
    while(high - low > 1) {
        int mid = low + ((high - low) / 2);
        if(i < mid) {
            // go left
            high = mid;
            
            if(cur_node->left == nullptr) {
                cur_node->left = new KPVector::Node(cur_node);
                cur_node = cur_node->left;
            }
        }
        else {
            // go right
            low = mid;

            if(cur_node->right == nullptr) {
                cur_node->right = new KPVector::Node(cur_node);
                cur_node = cur_node->right;
            }
        }
    }

    // At this point high - low == 1, so i == low
    assert(i == low);
    return cur_node;
}

double KPVector::get(int i) {
    assert(0 <= i && i < dim);

    Node* leaf = find(i);
    assert(leaf->sign != 0); // make sure it's actually a leaf

    double value = std::sqrt(leaf->weight);    
    if(leaf->sign == -1) {
        value = -value;
    }
    return value;
}

void KPVector::set(int i, double v) {
    assert(0 <= i && i < dim);

    Node* leaf = find(i);
    double old_weight = leaf->weight;
    leaf->weight = v * v;
    leaf->sign = (v < 0) ? -1 : 1;

    // Update parent weights 
    for(Node* cur_node = leaf->parent; cur_node != nullptr; cur_node = cur_node->parent) {
        cur_node->weight += (leaf->weight - old_weight);        
    }
}

double KPVector::squared_norm() {
    return root->weight;
}

double KPVector::sample() {
    Node* cur_node = root;
    std::uniform_real_distribution<double> unif_dist(0, 1);
    std::default_random_engine rand_eng;

    while(cur_node->sign == 0) {
        if(cur_node->left == nullptr && cur_node->right == nullptr) {
            throw std::logic_error("Vector is empty -- cannot sample.");
        }
        else if(cur_node->weight == 0) {
            return 0;
        }
        else if(cur_node->left == nullptr) {
            cur_node = cur_node->right;
        } 
        else if(cur_node->right == nullptr) {
            cur_node = cur_node->left;
        }
        else {
            double rand = unif_dist(rand_eng);
            if(rand <= (cur_node->left->weight / cur_node->weight)) {
                cur_node = cur_node->left;
            }
            else {
                cur_node = cur_node->right;
            }
        }
    }

    double value = std::sqrt(cur_node->weight);    
    if(cur_node->sign == -1) {
        value = -value;
    }
    return value;
}
