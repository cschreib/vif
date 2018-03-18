#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string sra = argv[1];
    std::string sdec = argv[2];
    double ra, dec;
    if (sex2deg(sra, sdec, ra, dec)) {
        print(format::precision(ra, 12), "   ", format::precision(dec, 12));
    } else {
        error("could not read sexagesimal coordinates");
        return 1;
    }

    return 0;
}

void print_help() {
    print("sex2deg v1.0");
    print("usage: sex2deg RA Dec");
}
