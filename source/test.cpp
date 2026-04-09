#include <iostream>
#include <string>
int main() {
    std::string g_algorithm = "todd_exp_000";
    std::cout << g_algorithm.length() << std::endl;
    std::cout << (g_algorithm.length() == 12) << std::endl;
    std::cout << g_algorithm.substr(0, 9) << std::endl;
    std::cout << (g_algorithm.substr(0, 9) == "todd_exp_") << std::endl;
    return 0;
}
