#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_help();
        return 0;
    }

    std::string ra = "";
    std::string dec = "";
    vec1s show;
    uint_t nsrc = 10;

    read_args(argc-3, argv+3, arg_list(ra, dec, nsrc, show));

    if (ra.empty() != dec.empty()) {
        error("you only specified 'ra' or 'dec', please provide both");
        return 1;
    }

    vec1d cra, cdec;

    if (!ra.empty()) {
        fits::read_table(argv[1], ra, cra, dec, cdec);
    } else {
        fits::read_table_loose(argv[1], "ra", cra, "dec", cdec);
        if (cra.empty()) {
            fits::read_table_loose(argv[1], "pos.ra", cra, "pos.dec", cdec);
            if (cra.empty()) {
                error("no standard position variable in this catalog (RA/DEC or POS.RA/POS.DEC)");
                note("please specify the correct variables using 'ra=...' and 'dec=...'");
                return 1;
            }
        }
    }

    vec2d show_data(show.size(), cra.size());

    for (uint_t i = 0; i < show.size(); ++i) {
        vec1d data;
        fits::read_table_loose(argv[1], show[i], data);
        if (data.empty()) {
            warning("no column named '", show[i], "', skipping");
            show[i] = "";
        }

        show_data(i,_) = data;
    }

    print("Catalog ranges:");
    print(" RA: [", min(cra), ", ", max(cra), "]");
    print(" Dec: [", min(cdec), ", ", max(cdec), "]");

    vec1d tra(1);
    vec1d tdec(1);
    std::string sra = argv[2];
    std::string sdec = argv[3];
    if (find(sra, ":") != npos) {
        note("using sexagesimal coordinates (hh:mm:ss, deg:mm:ss)");
        sex2deg(sra, sdec, tra[0], tdec[0]);
    } else {
        note("using decimal coordinates (degree)");
        if (!from_string(sra, tra[0])) {
            error("could not parse RA coordinate ", sra);
            return 1;
        }
        if (!from_string(sdec, tdec[0])) {
            error("could not parse Dec coordinate ", sdec);
            return 1;
        }
    }

    auto res = qxmatch(tra, tdec, cra, cdec, keywords(_nth(nsrc)));

    for (uint_t i = 0; i < nsrc; ++i) {
        std::string data;
        for (uint_t j = 0; j < show.size(); ++j) {
            if (show[j].empty()) continue;
            data += ", "+show[j]+"="+strn(show_data(j,res.id(i,0)));
        }
        print(i, ": ", res.d(i,0), "\", id=", res.id(i,0), data);
    }

    return 0;
}

void print_help() {
    print("find_source v1.0");
    print("usage: find_source catalog.fits 52.1456 -27.5465 [options]");
    print("");
    print("Available options:");
    print(" - ra: name of the Right Ascencion variable in the catalog");
    print(" - dec: name of the Declination variable in the catalog");
    print(" - show: array of column names to display");
    print(" - nsrc: number of nearby sources to show");
}
