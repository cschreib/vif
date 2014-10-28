#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    if (argc <= 2) {
        print_help();
        return 0;
    }

    struct cat_t {
        vec1d ra, dec;
    };

    cat_t cat;
    fits::read_table(argv[1], cat);

    cat_t fcat;
    fits::read_table(argv[2], fcat);

    vec1d range = {1.0, 40.0};
    uint_t nbin = 10;
    std::string out_file = "angcorrel.fits";
    uint_t tseed = 42;

    read_args(argc-2, argv+2, arg_list(range, nbin, name(out_file, "out"), name(tseed, "seed")));

    vec2d bins = e10(make_bins(log10(range[0]), log10(range[1]), nbin));
    vec1d ang = 0.5*(bins(0,_) + bins(1,_));

    vec1u hull = convex_hull(fcat.ra, fcat.dec);

    auto in_hull = [&](const vec1d& ra, const vec1d& dec) {
        return in_convex_hull(ra, dec, hull, fcat.ra, fcat.dec);
    };

    vec1d rra, rdec;
    randpos_uniform_options options;
    options.nsrc = cat.ra.size();
    if (options.nsrc < 1000) options.nsrc = 1000;

    auto seed = make_seed(tseed);
    auto status = randpos_uniform(seed,
        {min(fcat.ra[hull]),  max(fcat.ra[hull])},
        {min(fcat.dec[hull]), max(fcat.dec[hull])},
        in_hull, rra, rdec, options
    );

    if (!status.success) {
        error("could not generate random positions");
        note(status.failure);
        return 1;
    }

    vec1d w = angcorrel(cat.ra, cat.dec, rra, rdec, bins);

    fits::write_table(out_file, ftable(bins, w, ang));

    return 0;
}

void print_help() {
    using namespace format;

    print("angcorrel v1.0");
    paragraph("usage: angcorrel cat.fits refcat.fits [range,nbin,out,seed]");
    paragraph("Compute the angular two point correlation function of a given catalog "
        "'cat.fits'. The correlation is calculated using the Landy-Szalay estimator, by "
        "comparing against a random uniform distribution of points generated within the "
        "boundaries of the catalog 'refcat.fits'. The correlation function is written in "
        "a FITS file as column 'W', in units of counts per arcsec.");

    header("List of available command line options:");
    bullet("range", "[float,float] angular range within which to compute the correlation "
        "function [arcsec] (default: [1,40])");
    bullet("nbin", "[unsigned integer] number of angular bins (default: 10)");
    bullet("out", "[string] output file name (default: angcorrel.fits)");
    bullet("seed", "[unsigned integer] random seed for the generation of the random "
        "uniform positions (default: 42)");
}
