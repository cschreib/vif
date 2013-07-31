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
    vec1i gidg = where(finite(cat.flux) && finite(cat.flux_err) && cat.flux > 0 && cat.flux_err > 0
        && cat.flux/cat.flux_err > 3);

    std::string hband = "band";
    std::string hnote = "note";
    std::string hlam  = "lambda [um]";
    std::string hdepuJy = "depth [uJy]";
    std::string hdepAB = "depth [AB]";
    std::string hdet = "detections";
    std::string hd3s = "3-sigma";
    std::string hd5s = "5-sigma";

    uint_t maxb = std::max(hband.size(), max(length(cat.bands)));
    uint_t maxn = std::max(hnote.size(), max(length(cat.notes)));
    uint_t maxl = std::max(hlam.size(), max(length(strna(cat.lambda))));
    uint_t maxdJ = std::max(hdepuJy.size(), max(length(strna_sci(cat.flux_err[gidg]))));
    uint_t maxdAB = std::max(hdepAB.size(), max(length(strna(cat.flux_err[gidg]))));
    uint_t maxdet = std::max(hdet.size(), length(strn(ngal)));
    uint_t maxd3 = std::max(hd3s.size(), length(strn(ngal)));
    uint_t maxd5 = std::max(hd5s.size(), length(strn(ngal)));

    std::string header = " "+
        align_center(hband, maxb)+" | "+
        align_center(hnote, maxn)+" | "+
        align_center(hlam, maxl)+" | "+
        align_center(hdepuJy, maxdJ)+" | "+
        align_center(hdepAB, maxdAB)+" | "+
        align_center(hdet, maxdet)+" | "+
        align_center(hd3s, maxd3)+" | "+
        align_center(hd5s, maxd5)+" ";

    print(header);
    print(std::string(header.size(), '='));

    for (uint_t i = 0; i < n_elements(cat.bands); ++i) {
        auto f = cat.flux(i,_);
        auto e = cat.flux_err(i,_);

        vec1i idg = where(finite(f) && finite(e) && f > 0 && e > 0);

        print(" ",
            align_center(cat.bands[i], maxb), " | ",
            align_center(cat.notes[i], maxn), " | ",
            align_right(strn(cat.lambda[i]), maxl), " | ",
            align_right(strn_sci(median(e[idg])), maxdJ), " | ",
            align_right(strn(uJy2mag(median(e[idg]))), maxdAB), " | ",
            align_right(strn(total(f > 0)), maxdet), " | ",
            align_right(strn(total(f/e > 3)), maxd3), " | ",
            align_right(strn(total(f/e > 5)), maxd5)
        );
    }

    return 0;
}
