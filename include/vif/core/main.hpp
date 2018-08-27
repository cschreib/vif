#ifndef VIF_CORE_MAIN_HPP
#define VIF_CORE_MAIN_HPP

#include "vif/core/error.hpp"

#ifndef NO_LIBUNWIND
namespace vif {
    const char* executable_path;
}
#endif

// New entry point
int vif_main(int argc, char* argv[]);

// Define the standard main function
int main(int argc, char* argv[]) {
    #ifndef NO_LIBUNWIND
    vif::executable_path = argv[0];
    #endif

    return vif_main(argc, argv);
}

#endif
