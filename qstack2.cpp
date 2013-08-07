#include <phypp.hpp>

void print_help() {

}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string cat = "";
    std::string img = "";
    std::string wht = "";
    std::string err = "";
    std::string pos = "";
    std::string out = "";
    uint_t hsize = npos;
    bool median = false;
    bool mean = false;
    bool verbose = false;
    read_args(argc, argv, arg_list(out, cat, img, wht, err, pos, hsize, median, mean, verbose));

    if (!median && !mean) {
        mean = true;
    } else if (mean && median) {
        median = false;
    }

    if (hsize == npos) {
        error("missing cutout size 'hsize'");
        print_help();
        return 1;
    }

    if (out.empty()) {
        error("missing output file name 'out'");
        print_help();
        return 1;
    }

    if (cat.empty()) {
        error("missing catalog file 'cat'");
        print_help();
        return 1;
    } else if (!file::exists(cat)) {
        error("cannot find '"+cat+"'");
        return 1;
    }

    if (img.empty()) {
        error("missing image file 'img'");
        print_help();
        return 1;
    } else if (!file::exists(img)) {
        error("cannot find '"+img+"'");
        return 1;
    }

    if (!wht.empty() && !file::exists(wht)) {
        error("cannot find '"+wht+"'");
        return 1;
    }

    if (!err.empty() && !file::exists(err)) {
        error("cannot find '"+err+"'");
        return 1;
    }

    struct {
        vec1d ra, dec;
    } fcat;

    std::string posh = pos.empty() ? "" : pos+".";

    fits::read_table(cat, toupper(posh+"ra"), fcat.ra, toupper(posh+"dec"), fcat.dec);

    vec3f cube(0, (2*hsize+1), (2*hsize+1));
    vec1u fids;
    qstack(fcat.ra, fcat.dec, img, hsize, cube, fids);

    vec2f stack;
    if ((wht.empty() && err.empty()) || median) {
        if (verbose) print("stacking ", cube.dims[0], "/", fcat.ra.size(), " sources");
        if (mean) {
            stack = qstack_mean(cube);
        } else {
            stack = qstack_median(cube);
        }
    } else {
        if (median) {
            warning("cannot use weights when doing median stacking");
        }

        if (!wht.empty()) {
            vec3f wcube;
            vec1u wids;
            qstack(fcat.ra, fcat.dec, wht, hsize, wcube, wids);

            vec1u ids1, ids2;
            match(fids, wids, ids1, ids2);

            if (verbose) print("stacking ", ids1.size(), "/", fcat.ra.size(), " sources");

            stack = qstack_mean(cube(ids1,_,_), wcube(ids2,_,_));
        } else {
            vec3f ecube;
            vec1u eids;
            qstack(fcat.ra, fcat.dec, err, hsize, ecube, eids);

            vec1u ids1, ids2;
            match(fids, eids, ids1, ids2);

            if (verbose) print("stacking ", ids1.size(), "/", fcat.ra.size(), " sources");

            stack = qstack_mean(cube(ids1,_,_), pow(ecube(ids2,_,_), -2));
        }
    }

    file::mkdir(file::get_directory(out));
    fits::write(out, stack);

    return 0;
}
