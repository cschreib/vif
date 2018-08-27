#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    double ra, dec;
    if (!from_string(argv[1], ra)) {
        error("could not read RA");
        return 1;
    }
    if (!from_string(argv[2], dec)) {
        error("could not read Dec");
        return 1;
    }

    std::string sra, sdec;
    deg2sex(ra, dec, sra, sdec);
    print(sra, "   ", sdec);

    return 0;
}

void print_help() {
    print("deg2sex v1.0");
    print("usage: deg2sex RA Dec");
}
