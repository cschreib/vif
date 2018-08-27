#ifndef PHYPP_CORE_MAIN_HPP
#define PHYPP_CORE_MAIN_HPP

#include "phypp/core/error.hpp"

#ifndef NO_LIBUNWIND
namespace phypp {
    const char* executable_path;
}
#endif

// New entry point
int phypp_main(int argc, char* argv[]);

// Define the standard main function
int main(int argc, char* argv[]) {
    #ifndef NO_LIBUNWIND
    phypp::executable_path = argv[0];
    #endif

    return phypp_main(argc, argv);
}

#endif
