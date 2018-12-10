#include <iostream>
#include <fstream>

#include "KPMatrix.hh"
#include "ModFKV.hh"
#include "TangsAlgorithm.hh"

int main() {

    std::cout << "Program start" << std::endl;

    int m = 3;
    int n = 5;
    KPMatrix A(m, n);

    A.set(0, 1, 0.5);
    A.set(0, 2, 0.75);
    A.set(1, 0, 0.33);
    A.set(1, 1, 0.75);
    A.set(1, 2, 0.23);
    A.set(1, 3, 0.2);
    A.set(1, 4, 0.1);
    A.set(2, 0, 0.3);
    A.set(2, 1, 0.9);
    A.set(2, 2, 0.1);
    A.set(2, 3, 0.2);
    A.set(2, 4, 0.2);

    double threshold = 0.1; // between 0 and 1
    double epsilon = .01;
    double kappa = 0.5; // between 0 and 1
    TangSampler ts(A, threshold, epsilon, kappa);
    ts.sample(0);

    /*
    int m = 63978;
    int n = 150;
    KPMatrix A(m, n);

    // Read in data
    std::ifstream data_file("jester_ratings.dat");
    assert(data_file.is_open());

    int user;
    int item;
    double value;
    while(data_file.eof() == false) {
        data_file >> user;
        data_file >> item;
        data_file >> value;

        // Zero-indexed, not one-indexed
        A.set(user - 1, item - 1, value);
    }

    double threshold = 0.5; // between 0 and 1
    double epsilon = .001;
    double kappa = 0.5; // between 0 and 1
    TangSampler ts(A, threshold, epsilon, kappa);
    // Get a recommendation for user 0
    ts.sample(0);
    */
}
