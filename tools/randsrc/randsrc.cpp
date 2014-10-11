#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string pos = "";
    std::string out;
    uint_t tseed = 42;
    uint_t max_iter = 1000;
    uint_t nsrc = 0;
    read_args(argc-1, argv+1, arg_list(pos, out, nsrc, max_iter, name(tseed, "seed")));

    auto seed = make_seed(tseed);

    vec1d hra, hdec;
    if (!pos.empty() && !end_with(pos, ".")) pos = pos+".";
    fits::read_table(argv[1], pos+"ra", hra, pos+"dec", hdec);

    if (nsrc == 0) nsrc = hra.size();

    if (hra.size() < 3) {
        error("too few points to build convex hull");
        note("need at least 3");
        return 1;
    }

    vec1u hull = convex_hull(hra, hdec);

    vec1d rra = {min(hra), max(hra)};
    vec1d rdec = {min(hdec), max(hdec)};

    vec1d ra = randomu(seed, nsrc)*(rra[1] - rra[0]) + rra[0];
    vec1d dec = randomu(seed, nsrc)*(rdec[1] - rdec[0]) + rdec[0];

    vec1u bid = where(!in_convex_hull(ra, dec, hull, hra, hdec));
    uint_t iter = 0;
    while (!bid.empty() && iter < max_iter) {
        ra[bid] = randomu(seed, bid.size())*(rra[1] - rra[0]) + rra[0];
        dec[bid] = randomu(seed, bid.size())*(rdec[1] - rdec[0]) + rdec[0];
        bid = bid[where(!in_convex_hull(ra[bid], dec[bid], hull, hra, hdec))];
        ++iter;
    }

    if (iter == max_iter) {
        warning("maximum number of iteration reached (", max_iter, ")");
        note("some sources will be out of the convex hull");
        note("re-run the program with a higher value of 'max_iter' to fix this");
    }

    if (out.empty()) {
        out = argv[1];
        if (end_with(out, ".fits")) {
            out = erase_end(out, ".fits");
        }

        out += "-rnd.fits";
    }

    file::mkdir(file::get_directory(out));

    fits::write_table(out, ftable(ra, dec));

    return 0;
}

void print_help() {
    print("randsrc v1.0");
    print("usage: randsrc cat.fits [seed,out,nsrc,pos,max_iter]");
}
