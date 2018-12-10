#pragma once

#include <random>

#include "KPMatrix.hh"
#include "ModFKV.hh"

int sample_from_distribution(const KPVector& vector, const ModFKV& rank_reduction) {
    assert(vector.dimension() == rank_reduction.get_rank());

    int k = vector.dimension();
    int n = rank_reduction.get_V_hat_column(0).dimension();

    KPVector distribution_P(k);
    for(int i = 0; i < k; ++i) {
        distribution_P.set(i, vector.get(i) * std::sqrt(rank_reduction.get_V_hat_column(i).squared_norm()));
    }

    // Establish random number generator
    std::uniform_real_distribution<double> unif_dist(0, 1);
    std::random_device rd;
    std::default_random_engine rand_eng(rd());

    // Sample until we don't reject one
    while(true) {
        // col is in [k]
        int col = distribution_P.sample_index();
        // sample is in [n]
        int sample = rank_reduction.get_V_hat_column(col).sample_index();

        // Compute probability of accepting sample
        double denominator = 0;
        for(int j = 0; j < k; ++j) {
            double wV_j = vector.get(j) * rank_reduction.get_V_hat_column(j).get(sample);
            denominator += wV_j * wV_j;
        }
        denominator *= vector.dimension();

        double numerator = 0;
        // Calculate (Vw)^2_sample
        for(int j = 0; j < k; ++j) {
            numerator += vector.get(j) * rank_reduction.get_V_hat_column(j).get(sample);
        }
        numerator = std::abs(numerator);
        denominator = std::abs(denominator);

        // assert(0 < numerator && numerator <= denominator);
        std::cout << "Probability of acceptance of sample: " << numerator << " / " << denominator << std::endl;

        // Accept or reject sample     
        if(unif_dist(rand_eng) < numerator / denominator) {
            return sample;
        }
    }
}

