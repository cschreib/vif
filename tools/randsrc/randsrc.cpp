#include <phypp.hpp>

void print_help();

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

    // Compute bounding box of input positions
    vec1d rra = {min(hra), max(hra)};
    vec1d rdec = {min(hdec), max(hdec)};

    vec1s vmodel = split(model, ":");

    vec1d ra, dec;

    if (vmodel[0] == "uniform") {
        // Place points uniformly within the convex hull formed by the provided coordinates
        ra = (randomu(seed, nsrc) - 0.5)*(rra[1] - rra[0]) + mean(rra);
        dec = (randomu(seed, nsrc) - 0.5)*(rdec[1] - rdec[0]) + mean(rdec);

        vec1u bid = where(!in_convex_hull(ra, dec, hull, hra, hdec));
        uint_t iter = 0;
        while (!bid.empty() && iter < max_iter) {
            ra[bid] = (randomu(seed, bid.size()) - 0.5)*(rra[1] - rra[0]) + mean(rra);
            dec[bid] = (randomu(seed, bid.size()) - 0.5)*(rdec[1] - rdec[0]) + mean(rdec);
            bid = bid[where(!in_convex_hull(ra[bid], dec[bid], hull, hra, hdec))];
            ++iter;
        }

        if (iter == max_iter) {
            warning("maximum number of iteration reached (", max_iter, ")");
            note("some sources will be out of the convex hull");
            note("re-run the program with a higher value of 'max_iter' to fix this");
        }
    } else if (vmodel[0] == "power") {
        // Place points with a power law angular correlation function within the convex
        // hull formed by the provided coordinates, using the Soneira & Peebles algorithm.
        // http://www.astro.rug.nl/~weygaert/tim1publication/lss2007/soneirapeebles.pdf

        // Read the power law index in 'omega(theta) ~ r^(-power)'
        double power;
        if (vmodel.size() <= 1) {
            power = 1.0;
        } else if (!from_string(vmodel[1], power)) {
            error("could not read power law index from '", vmodel[1], "'");
            return 1;
        }

        // Read the number of levels to use in the generator
        uint_t nlevel;
        if (vmodel.size() <= 2) {
            nlevel = 5;
        } else if (!from_string(vmodel[2], nlevel)) {
            error("could not read number of levels from '", vmodel[2], "'");
            return 1;
        }

        // Choose the number of objects to generate
        // If the wanted number of object is too low, the algorithm will not work
        // well, so we generate more and will randomly pick among those afterwards.
        uint_t nsim = std::max(nsrc, uint_t(1000));

        // Compute initial parameters for the algorithm
        // The initial radius is taken from the diagonals of the RA/Dec bounding box.
        double r0 = std::max(
            angdist(rra[0], rdec[0], rra[1], rdec[1]),
            angdist(rra[0], rdec[1], rra[1], rdec[0])
        )/3600.0/2.0;

        // The number of object per cell
        uint_t eta = ceil(pow(nsim, 1.0/nlevel));
        if (eta <= 1) {
            eta = 2;
            nlevel = ceil(log10(nsim)/log10(eta));
        }

        // The radius shrinking factor
        double lambda = pow(eta, 1.0/(2.0 - power));

        if (verbose) {
            print("power: ", power, ", nlevel: ", nlevel,
                ", eta: ", eta, ", lambda: ", lambda, ", r0=", r0);
        }

        auto get_fill = [&hra, &hdec, &hull, &seed, max_iter](double ra0, double dec0,
            double r, uint_t num) mutable {

            vec1d tdec = (randomu(seed, num) - 0.5)*2*r + dec0;
            vec1d tra  = (randomu(seed, num) - 0.5)*2*r/cos(tdec*dpi/180.0) + ra0;

            vec1u id = where(angdist_less(tra, tdec, ra0, dec0, r*3600.0));

            return fraction_of(in_convex_hull(tra[id], tdec[id], hull, hra, hdec));
        };

        auto gen = [&hra, &hdec, &hull, &seed, max_iter](double ra0, double dec0,
            double r, uint_t num, vec1d& tra, vec1d& tdec) mutable {

            tdec = (randomu(seed, num) - 0.5)*2*r + dec0;
            tra  = (randomu(seed, num) - 0.5)*2*r/cos(tdec*dpi/180.0) + ra0;

            vec1u bid = where(
                !in_convex_hull(tra, tdec, hull, hra, hdec) ||
                !angdist_less(tra, tdec, ra0, dec0, r*3600.0)
            );

            uint_t iter = 0;
            while (!bid.empty() && iter < max_iter) {
                tdec[bid] = (randomu(seed, bid.size()) - 0.5)*2*r + dec0;
                tra[bid]  = (randomu(seed, bid.size()) - 0.5)*2*r/cos(tdec[bid]*dpi/180.0) + ra0;
                bid = bid[where(
                    !in_convex_hull(tra[bid], tdec[bid], hull, hra, hdec) ||
                    !angdist_less(tra[bid], tdec[bid], ra0, dec0, r*3600.0)
                )];
                ++iter;
            }

            if (iter == max_iter) {
                warning("maximum number of iteration reached (", max_iter, ")");
                note("some sources will be out of the convex hull");
                note("re-run the program with a higher value of 'max_iter' to fix this");
            }
        };

        // Initialize the algorithm with random positions in the whole area
        gen(mean(rra), mean(rdec), r0, eta, ra, dec);

        for (int_t l : range(nlevel-1)) {
            // For each position, generate new objects within a "cell" of radius 'r1'
            double r1 = r0*pow(lambda, -l-1);

            // First, compute the filling factor of each cell, i.e. what fraction of the
            // area of the cell is inside the convex hull. Each cell will then only
            // generate this fraction of object, in order not to introduce higher densities
            // close to the borders of the field.
            vec1d fill(ra.size());
            for (uint_t i : range(ra)) {
                fill[i] = get_fill(ra[i], dec[i], r1, 100u);
            }

            // Normalize filling factors so that the total number of generated object
            // is kept constant.
            fill /= mean(fill);

            // Assign a number of objects to each positions
            vec1u ngal = floor(eta*fill);
            vec1d dngal = eta*fill - ngal;
            uint_t miss = pow(eta, l+2) - total(ngal);

            // Randomly assign fractional number of objects to make sure that the total
            // number of object is preserved.
            uint_t iter = 0;
            while (miss != 0 && iter < max_iter) {
                for (uint_t i : range(ra)) {
                    if (dngal[i] != 0.0 && randomu(seed) < dngal[i]) {
                        ++ngal[i];
                        dngal[i] = 0.0;

                        --miss;
                        if (miss == 0) break;
                    }
                }

                ++iter;
            }

            // Generate new positions
            vec1d ra1, dec1;
            for (uint_t i : range(ra)) {
                if (ngal[i] == 0) continue;

                vec1d tra1, tdec1;
                gen(ra[i], dec[i], r1, ngal[i], tra1, tdec1);

                append(ra1, tra1);
                append(dec1, tdec1);
            }

            // Use the newly generated positions as input for the next level
            std::swap(ra, ra1);
            std::swap(dec, dec1);
        }

        // Trim catalog to match the input number of sources
        vec1u fids = shuffle(uindgen(ra.size()), seed)[uindgen(nsrc)];
        ra = ra[fids];
        dec = dec[fids];
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
