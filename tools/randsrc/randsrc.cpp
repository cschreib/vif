#include <phypp.hpp>

void print_help();

template<typename T, typename D, typename I>
bool read_param(T& p, D def, std::string name, I& iter, const I& end) {
    if (iter == end) {
        p = def;
    } else if (!from_string(*iter, p)) {
        error("could not read ", name, " from '", *iter, "'");
        return false;
    } else {
        ++iter;
    }

    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    // Program arguments
    std::string pos = "";
    std::string out;
    uint_t tseed = 42;
    uint_t max_iter = 1000;
    uint_t nsrc = 0;
    std::string model = "uniform";
    bool verbose = false;

    read_args(argc-1, argv+1, arg_list(pos, out, nsrc, max_iter, model,
        name(tseed, "seed"), verbose));

    auto seed = make_seed(tseed);

    // Read input position list
    vec1d hra, hdec;
    if (!pos.empty() && !end_with(pos, ".")) pos = pos+".";
    fits::read_table(argv[1], pos+"ra", hra, pos+"dec", hdec);

    if (nsrc == 0) nsrc = hra.size();

    if (hra.size() < 3) {
        error("too few points to build convex hull");
        note("need at least 3");
        return 1;
    }

    // Compute convex hull of input positions
    vec1u hull = convex_hull(hra, hdec);
    auto in_hull = [&](const vec1d& tra, const vec1d& tdec) {
        return in_convex_hull(tra, tdec, hull, hra, hdec);
    };

    // Compute bounding box of input positions
    vec1d rra = {min(hra), max(hra)};
    vec1d rdec = {min(hdec), max(hdec)};

    vec1s vmodel = split(model, ":");

    vec1d ra, dec;

    if (vmodel[0] == "uniform") {
        // Place points uniformly within the convex hull formed by the provided coordinates
        randpos_uniform_options opt;
        opt.nsrc = nsrc;
        opt.max_iter = max_iter;

        auto status = randpos_uniform(seed, rra, rdec, in_hull, ra, dec, opt);
        if (!status.success) {
            error("generation failed: ", status.failure);
            return 1;
        }
    } else if (vmodel[0] == "power") {
        // Place points with a power law angular correlation function within the convex
        // hull formed by the provided coordinates, using the Soneira & Peebles algorithm.
        // http://www.astro.rug.nl/~weygaert/tim1publication/lss2007/soneirapeebles.pdf

        randpos_power_options opt;
        opt.nsrc = nsrc;
        opt.max_iter = max_iter;

        auto iter = vmodel.begin()+1;

        // Read the power law index in 'omega(theta) ~ r^(-power)'
        if (!read_param(opt.power, 1.0, "power law index", iter, vmodel.end())) {
            return 1;
        }

        // Read the number of levels to use in the generator
        if (!read_param(opt.levels, 5, "number of levels", iter, vmodel.end())) {
            return 1;
        }

        // Read the number of sources to generate and pick from
        if (!read_param(opt.overload, 1.0, "overload factor", iter, vmodel.end())) {
            return 1;
        }

        // Read the area inflation parameter
        if (!read_param(opt.inflate, 1.0, "inflation factor", iter, vmodel.end())) {
            return 1;
        }

        if (vmodel.size() > 5) {
            warning("too many parameters provided (maximum is 4)");
        }

        auto status = randpos_power(seed, hull, hra, hdec, ra, dec, opt);
        if (!status.success) {
            error("generation failed: ", status.failure);
            return 1;
        }
    } else {
        error("unknown model '", vmodel[0], "'");
        return 1;
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
    print("usage: randsrc cat.fits [seed,out,nsrc,pos,max_iter,model,verbose]");
}
