#pragma once

#include "Eigen/Dense"
#include "KPMatrix.hh"
#include "ModFKV.hh"


class TangSampler {
    public:
        TangSampler(const KPMatrix& A, double threshold, double epsilon, double kappa) 
            : A_(A)    
        {
            // Bounds check
            assert(0 < epsilon);
            assert(0 < kappa && kappa <= 1);
            assert(0 < threshold);

            rank_reduction_ = ModFKV(A, 
                    threshold * (1 - (kappa / 2)),
                    epsilon,
                    2 * kappa / (2 - kappa));
        }

        // Sample for a given user
        int sample(int user) {
            assert(0 <= user && user < A_.num_rows());
            // TODO: implement this
            assert(false);
        }

    private:
        const KPMatrix& A_;
        ModFKV rank_reduction_;
        

};
