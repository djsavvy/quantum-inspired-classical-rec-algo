#pragma once

#include "Eigen/Dense"
#include "KPMatrix.hh"
#include "ModFKV.hh"
#include "SampleFromDist.hh"


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
            return sample_from_distribution(est, rank_reduction_);
        }

    private:
        const KPMatrix& A_;
        ModFKV rank_reduction_;
        double epsilon_;
        double kappa_;

};
