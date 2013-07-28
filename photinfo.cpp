#include <phypp.hpp>
#include <filters.hpp>

int main(int argc, char* argv[]) {
    if (argc == 1) {
        return 0;
    }

    read_args(argc-1, argv+1, arg_list(data_dir));

    struct {
        vec2f flux, flux_err;
        vec1s bands, notes;
        vec1f lambda;
    } cat;

    fits::read_table_loose(argv[1], cat);

    uint_t ngal = dim(cat.flux)[1];

    if (n_elements(cat.lambda) == 0) {
        cat.lambda = fltarr(n_elements(cat.bands));
        for (uint_t i = 0; i < n_elements(cat.bands); ++i) {
            auto filter = get_filter(cat.bands[i]);
            cat.lambda[i] = filter.rlam;
        }
    }

    if (n_elements(cat.notes) == 0) {
        cat.notes = replicate("", n_elements(cat.bands));
    }

    print("Photometry: [", ngal, " sources]");
    for (uint_t i = 0; i < n_elements(cat.bands); ++i) {
        print("- ", cat.bands[i], " (", cat.notes[i], "): ", cat.lambda[i], "um");
        print("  ", strn(total(cat.flux(i,_) > 0), floor(log10(ngal))+1, ' '), " detections",
            ", ", strn(total(cat.flux(i,_)/cat.flux_err(i,_) > 3), floor(log10(ngal))+1, ' '),
            " (3-sigma)",
            ", ", strn(total(cat.flux(i,_)/cat.flux_err(i,_) > 5), floor(log10(ngal))+1, ' '),
            " (5-sigma)"
        );
    }

    return 0;
}
