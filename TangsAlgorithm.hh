#pragma once

#include "Eigen/Dense"
#include "KPMatrix.hh"
#include "ModFKV.hh"


class TangSampler {
    public:
        TangSampler(const KPMatrix& A, double threshold, double epsilon, double kappa) 
            : A_(A)
            , rank_reduction_(A,
                    threshold * (1 - kappa / 2),
                    epsilon,
                    2 * kappa / (2 - kappa)
              ) 
            , epsilon_(epsilon)
            , kappa_(kappa)
        {
            // Bounds check
            assert(0 < epsilon);
            assert(0 < kappa && kappa <= 1);
            assert(0 < threshold);
        }

        // Sample for a given user
        int sample(int user) {
            assert(0 <= user && user < A_.num_rows());

            // Calculate est
            int k = rank_reduction_.get_rank();
            KPVector est(k);
            for(int t = 0; t < k; ++t) {
                est.set(k, estimate_dot_product(
                            A_.get_row(user), 
                            rank_reduction_.get_V_hat_column(t),
                            epsilon_,
                            epsilon_ / std::sqrt(k)
                        ));
            }

            // Sample from est V_hat^T



        }

    private:
        const KPMatrix& A_;
        ModFKV rank_reduction_;
        double epsilon_;
        double kappa_;
        

};

int sample_from_distribution(const KPVector& vector, const ModFKV& rank_reduction) {
    assert(vector.dimension() == rank_reduction.get_rank());

    KPVector distribution_P(vector.dimension());
    for(int i = 0; i < vector.dimension(); ++i) {
        distribution_P.set(i, vector.get(i) * std::sqrt(rank_reduction.get_V_hat_column(i).squared_norm()));
    }

    // Sample until we don't reject one
    while(true) {
        int col = distribution_P.sample_index();
        int sample = rank_reduction.get_V_hat_column(col).sample_index();

        double denominator = 0;
        for(int j = 0; j < vector.dimension(); ++j) {
            double wV_j = vector.get(j) * rank_reduction.get_V_hat_column(j).get(col);
            denominator += wV_j * wV_j;
        }
        denominator *= vector.dimension();

        double numerator = 0;
        for(int j = 0; j < ; ++j) {
            numerator += vector.get(j) * rank_reduction.get_V_hat_column(col).get(j);
        }

        
        

    }
}
