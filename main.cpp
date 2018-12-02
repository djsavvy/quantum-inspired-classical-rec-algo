#include <iostream>
#include "KPMatrix.hh"

int main() {

    KPMatrix test_matr(5, 5);
    test_matr.set(0, 2, -3);
    test_matr.set(3, 4, 5);
    test_matr.set(4, 0, -1);

    std::cout << test_matr.get(0, 2) << std::endl;
    std::cout << test_matr.get(3, 4) << std::endl;
    std::cout << test_matr.squared_frobenius_norm() << std::endl;
    std::cout << test_matr.sample_a_row() << std::endl;
    std::cout << test_matr.sample_from_row(0) << std::endl;
    std::cout << test_matr.sample_from_row(2) << std::endl;



    /*
    KPVector test_vec(10);
    test_vec.set(0, 2);
    test_vec.set(6, 7);
    test_vec.set(1, -3);

    std::cout << test_vec.get(0) << std::endl;
    std::cout << test_vec.get(1) << std::endl;
    std::cout << test_vec.get(6) << std::endl;

    std::cout << "----" << std::endl;

    std::cout << test_vec.squared_norm() << std::endl;

    std::cout << "----" << std::endl;

    for(int i = 0; i < 10; ++i) {
        std::cout << "Sample value: " << test_vec.get(test_vec.sample_index()) << std::endl;
        std::cout << "Sample index: " << test_vec.sample_index() << std::endl;
    }

    std::cout << "----" << std::endl;

    std::cout << test_vec << std::endl;
    */

}
