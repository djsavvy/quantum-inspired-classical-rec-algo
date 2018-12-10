#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
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
            if(p <= 0 || p > A.num_rows()) {
                p = A.num_rows();
                // std::cout << p << std::endl;
            }
            
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

            // std::cout << "Computed W" << std::endl;

            // Compute Singular Value Decomposition of W
            Eigen::BDCSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeThinU);
            new_rank = 0;
            while(new_rank < svd.singularValues().size() 
                    && svd.singularValues()(new_rank) > threshold) {
                ++new_rank;
            }
            if(new_rank <= 0) new_rank = 1;

            singular_values = svd.singularValues();
            singular_values.resize(p);
            left_singular_vectors = svd.matrixU();
            left_singular_vectors.resize(Eigen::NoChange, p);

            // std::cout << "Done with SVD" << std::endl;

            // Compute columns of the V_hat representation for use in the main algorithm
            V_hat_columns.reserve(p);
            for(int i = 0; i < p; ++i) {
                // std::cout << "Starting i loop " << i << std::endl;
                V_hat_columns.push_back(KPVector(A.num_cols()));
                for(int j = 0; j < A.num_cols(); ++j) {
                    double value = 0;
                    for(int x = 0; x < p; ++x) {
                        value += A.get(row_indices[x], j) * left_singular_vectors(x, i);
                        // std::cout << "got value " << x << " " << j << " " << i << std::endl;
                    }
                    value /= singular_values(i);
                    // std::cout << "got singular value " << i << std::endl;

                    V_hat_columns[i].set(j, value);
                    // std::cout << "done with j loop " << j << std::endl;
                }
            }
            // std::cout << "Done constructing ModFKV" << std::endl;
        }

        int get_rank() const {
            // return p;
            // std::cout << "V_hat_columns.size(): " << V_hat_columns.size() << std::endl;
            return V_hat_columns.size();
        }

        const KPVector& get_V_hat_column(int index) const {
            assert(0 <= index && index < V_hat_columns.size());
            return V_hat_columns[index];
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

        std::vector<KPVector> V_hat_columns;
};
