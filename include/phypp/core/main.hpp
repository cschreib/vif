#ifndef PHYPP_CORE_MAIN_HPP
#define PHYPP_CORE_MAIN_HPP

#include "phypp/core/error.hpp"

const char* executable_path;

int phypp_main(int argc, char* argv[]);

// Define the standard main function
int main(int argc, char* argv[]) {
    executable_path = argv[0];
    return phypp_main(argc, argv);
}

#endif
