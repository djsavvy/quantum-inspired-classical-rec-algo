#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <set>

#include "KPMatrix.hh"
#include "Eigen/Dense"


class ModFKV {
    // Finds a low-rank approximation for a given KPMatrix
    
    public:
        ModFKV(const KPMatrix& A, double threshold, double epsilon, double kappa) 
        {
            // Bounds check
            double frobenius_norm = std::sqrt(A.squared_frobenius_norm());
            assert( 0 < threshold && threshold <= frobenius_norm);
            assert(
                0 < epsilon 
                && (16 * epsilon * epsilon) <= (threshold / frobenius_norm) 
            );
            assert( (epsilon * epsilon <= kappa) && (kappa <= 1) );

            // Compute constants
            K = A.squared_frobenius_norm() / (threshold * threshold);
            e_bar = kappa * epsilon * epsilon / std::sqrt(K);
            p = static_cast<int>(std::ceil(
                10000000 * std::max(
                    (K * K * K * K) / (e_bar * e_bar * e_bar),
                    (K * K) / (e_bar * e_bar * e_bar * e_bar)
                )
            ));
            assert( 0 < p && p <= A.num_rows() );
            
            // Sample rows from D_{A~}
            row_indices.reserve(p);
            for(int i = 0; i < p; ++i) {
                row_indices.push_back(A.sample_a_row());
            }

            // Compute the distribution F
            Eigen::VectorXd density_F(A.num_cols());
            density_F.setZero(A.num_cols());
            for(int i = 0; i < p; ++i) {
                for(int j = 0; j < A.num_cols(); ++j) {
                    double val_from_A = A.get(row_indices[i], j);
                    density_F(j) += (val_from_A * val_from_A / A.row_squared_norm(row_indices[i]));
                }
            }
            density_F /= p;

            // Sample columns
            col_indices.reserve(p);
            std::uniform_int_distribution<> unif_dist(0, p - 1);
            std::random_device rd;
            std::default_random_engine rand_eng(rd());
            for(int i = 0; i < p; ++i) {
                int row_to_sample_from = unif_dist(rand_eng);
                col_indices.push_back(A.sample_from_row(row_to_sample_from));
            }

            // Put sampled, scaled values into matrix W
            Eigen::MatrixXd W(p, p);
            for(int i = 0; i < p; ++i) {
                for(int j = 0; j < p; ++j) {
                    double scaled_value = A.get(row_indices[i], col_indices[j]) / p;
                    scaled_value /= std::sqrt( A.row_squared_norm(i) / A.squared_frobenius_norm() );
                    scaled_value /= std::sqrt( density_F(j) );
                    W(i, j) = scaled_value;
                }
            }

            // Compute Singular Value Decomposition of W
            Eigen::BDCSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeThinU);
            new_rank = 0;
            while(new_rank < svd.singularValues().size() 
                    && svd.singularValues()(new_rank) > threshold) {
                ++new_rank;
            }
            assert(new_rank > 0);

            singular_values = svd.singularValues();
            singular_values.resize(new_rank);
            left_singular_vectors = svd.matrixU();
            left_singular_vectors.resize(Eigen::NoChange, new_rank);
        }

    private:
        double K;
        double e_bar;
        int p;

        std::vector<int> row_indices;
        std::vector<int> col_indices;

        int new_rank;
        Eigen::VectorXd singular_values;
        Eigen::MatrixXd left_singular_vectors; // The vectors are the columns
};
