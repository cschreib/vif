#include <vif.hpp>
#include <vif/astro/qxmatch.hpp>

using namespace vif;
using namespace vif::astro;

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 4) {
        print_help();
        return 0;
    }

    std::string ra = "";
    std::string dec = "";
    vec1s show;
    uint_t nsrc = 10;
    std::string region = "";
    vec1s region_text;

    read_args(argc-3, argv+3, arg_list(ra, dec, nsrc, show, region, region_text));

    vec1u id; vec1d d;

    vec1d cra, cdec;

    if (ra.empty() != dec.empty()) {
        error("you only specified 'ra' or 'dec', please provide both");
        return 1;
    }

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

    vec1d tra(1);
    vec1d tdec(1);
    bool xmatch = true;

    if (std::string(argv[2]) == "id") {
        vec1s vs = split(trim(argv[3], "[]"), ",");
        vec1b parsed = from_string(vs, id);
        if (total(!parsed) != 0) {
            error("could not parse source IDs: ", ("'"+vs[where(!parsed)]+"'"));
            return 1;
        }

        if (id.size() > 1) {
            d.resize(id.size());
            xmatch = false;
        } else {
            tra[0] = cra[id[0]];
            tdec[0] = cdec[id[0]];
        }
    } else {
        print("Catalog ranges:");
        print(" RA: [", min(cra), ", ", max(cra), "]");
        print(" Dec: [", min(cdec), ", ", max(cdec), "]");

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
    }

    if (xmatch) {
        qxmatch_params p; p.nth = nsrc; p.brute_force = true;
        auto res = qxmatch(tra, tdec, cra, cdec, p);

        id = res.id(_,0);
        d = res.d(_,0);
    }

    vec2d show_data(show.size(), id.size());
    vec2d region_data(region_text.size(), id.size());

    for (uint_t i = 0; i < show.size(); ++i) {
        vec1d data;
        fits::read_table_loose(argv[1], show[i], data);
        if (data.empty()) {
            warning("no column named '", show[i], "', skipping");
            show[i] = "";
        } else {
            show_data(i,_) = data[id];
        }
    }

    for (uint_t i = 0; i < region_text.size(); ++i) {
        vec1d data;
        fits::read_table_loose(argv[1], region_text[i], data);
        if (data.empty()) {
            warning("no column named '", region_text[i], "', skipping");
            region_text[i] = "";
        } else {
            region_data(i,_) = data[id];
        }
    }

    for (uint_t i : range(id)) {
        std::string data;
        for (uint_t j = 0; j < show.size(); ++j) {
            if (show[j].empty()) continue;
            data += ", "+show[j]+"="+to_string(show_data(j,i));
        }
        print(i, ": ", d[i], "\", id=", id[i], data);
    }

    if (!region.empty()) {
        file::mkdir(file::get_directory(region));

        std::ofstream out(region);

        out << "# Region file format: DS9 version 4.1\n";
        out << "global color=green dashlist=8 3 width=2 font=\"helvetica 10 normal roman\""
            "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n";
        out << "fk5\n";

        for (uint_t i : range(id)) {
            std::string rra, rdec;
            deg2sex(cra[id[i]], cdec[id[i]], rra, rdec);
            out << "circle(" << rra << "," << rdec << ",0.5\") # width=3";
            if (!region_text.empty()) {
                out << " text={";
                bool first = true;
                for (uint_t j : range(region_text)) {
                    if (region_text[j].empty()) continue;
                    if (!first) out << ", ";
                    out << region_data(j,i);
                    first = false;
                }
                out << "}";
            } else {
                out << " text={" << id[i] << "}";
            }
            out << "\n";
        }
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
    print(" - region: output a DS9 region file with all the sources");
    print(" - region_text: value to display over each source in DS9");
}
