#include <iostream>
#include "KPMatrix.hh"

int main() {

    std::cout << "Hello world" << std::endl;

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
        std::cout << test_vec.sample() << std::endl;
    }

    std::cout << "----" << std::endl;

    std::cout << test_vec << std::endl;

}
