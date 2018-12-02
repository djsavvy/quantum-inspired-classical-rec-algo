#pragma once

#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>


class KPVector {
    public:
        KPVector(int n);
        ~KPVector();

        int dimension();
        double get(int i);
        void set(int i, double v);
        double squared_norm();
        int sample_index();

        friend std::ostream& operator<<(std::ostream& out, const KPVector& v);


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

            friend std::ostream& operator<<(std::ostream& out, const Node& n) {
                out << "Node: weight " << n.weight << ", sign " << n.sign << std::endl;
                if(n.left != nullptr) {
                    out << "Left: " << *(n.left) << std::endl;
                }
                if(n.right != nullptr) {
                    out << "Right: " << *(n.right) << std::endl;
                }
                return out;
            }
        };

        // Finds pointer to the leaf node for a given index
        // If it doesn't yet exist, creates it
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

int KPVector::dimension() {
    return dim;
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
        assert(high > low);
        int mid = low + ((high - low) / 2);
        if(i < mid) {
            // go left
            high = mid;
            
            if(cur_node->left == nullptr) {
                cur_node->left = new KPVector::Node(cur_node);
            }
            cur_node = cur_node->left;
        }
        else {
            // go right
            low = mid;

            if(cur_node->right == nullptr) {
                cur_node->right = new KPVector::Node(cur_node);
            }
            cur_node = cur_node->right;
        }
    }

    // At this point high - low == 1, so i == low
    assert(i == low);
    return cur_node;
}

double KPVector::get(int i) {
    assert(0 <= i && i < dim);

    Node* leaf = find(i);
    assert(leaf->sign != 0); // make sure it's actually a leaf, and has been set

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

int KPVector::sample_index() {
    Node* cur_node = root;
    std::uniform_real_distribution<double> unif_dist(0, 1);
    std::random_device rd;
    std::default_random_engine rand_eng(rd());

    int low = 0, high = dim;
    int mid = low + ((high - low) / 2);
    while(cur_node->sign == 0) {
        mid = low + ((high - low) / 2);
        if(cur_node->left == nullptr && cur_node->right == nullptr) {
            throw std::logic_error("Vector is empty -- cannot sample.");
        }
        else if(cur_node->left == nullptr) {
            cur_node = cur_node->right;
            low = mid;
        } 
        else if(cur_node->right == nullptr) {
            cur_node = cur_node->left;
            high = mid;
        }
        else {
            double rand = unif_dist(rand_eng);
            if(rand <= (cur_node->left->weight / cur_node->weight)) {
                cur_node = cur_node->left;
                high = mid;
            }
            else {
                cur_node = cur_node->right;
                low = mid;
            }
        }
    }

    int index = low;
    return index;
}

std::ostream& operator<<(std::ostream& out, const KPVector& v) {
    out << "Root: " << *(v.root);
    return out;
}



class KPMatrix {
    public:
        KPMatrix(int m, int n);
        ~KPMatrix();

        int num_rows();
        int num_cols();
        double get(int i, int j);
        void set(int i, int j, double v);
        double squared_frobenius_norm();
        // Get norm of a row
        double row_squared_norm(int i);

        // Sample over D_{A~}
        double sample_a_row();
        // Sample over D_{A_i}
        double sample_from_row(int i);

    private:
        int dim_m, dim_n;
        std::vector<KPVector*> vectors;
        double sum_vector_norms;

        // Norms of each row
        KPVector* vector_norms;
};


KPMatrix::KPMatrix(int m, int n) {
    assert(m > 0 && n > 0);
    dim_m = m;
    dim_n = n;

    vectors.reserve(m);
    for(int i = 0; i < m; ++i) {
        vectors.push_back(new KPVector(n));
    }
    
    vector_norms = new KPVector(n);
    sum_vector_norms = 0;
}

KPMatrix::~KPMatrix() {
    delete vector_norms;
    for(KPVector* v : vectors) {
        delete v;
    }
}

int KPMatrix::num_rows() {
    return dim_m;
}

int KPMatrix::num_cols() {
    return dim_n;
}

double KPMatrix::get(int i, int j) {
    assert(0 <= i && i < dim_m);
    assert(0 <= j && j < dim_n);
    return vectors[i]->get(j);
}

void KPMatrix::set(int i, int j, double v) {
    assert(0 <= i && i < dim_m);
    assert(0 <= j && j < dim_n);
    double old_row_norm = vectors[i]->squared_norm();
    vectors[i]->set(j, v);
    double new_row_norm = vectors[i]->squared_norm();
    vector_norms->set(i, new_row_norm);
    sum_vector_norms += (new_row_norm - old_row_norm);
}

double KPMatrix::squared_frobenius_norm() {
    return sum_vector_norms; 
}

double KPMatrix::row_squared_norm(int i) {
    assert(0 <= i && i < dim_m);
    return vectors[i]->squared_norm();
}

double KPMatrix::sample_a_row() {
    return vector_norms->sample_index();
}

double KPMatrix::sample_from_row(int i) {
    assert(0 <= i && i < dim_m);
    return vectors[i]->sample_index();
}
