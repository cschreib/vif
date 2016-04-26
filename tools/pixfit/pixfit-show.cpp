#include "pixfit-common.hpp"

int phypp_main(int argc, char* argv[]) {
    std::vector<map_info> maps;
    if (!read_maps(argv[3], maps)) return 1;

    std::string cat_file = argv[4];



    return 0;
}
