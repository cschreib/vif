#include <phypp.hpp>
#include <filters.hpp>

int main(int argc, char* argv[]) {
    if (argc == 1) {
        return 0;
    }

    read_args(argc-1, argv+1, arg_list(data_dir));

    struct {
        vec1d ra, dec;

        struct {
            vec1d ra, dec;
        } pos;

        vec2f flux, flux_err;
        vec1s bands, notes;
        vec1f lambda;
    } cat;

    fits::read_table_loose(argv[1], cat);

    uint_t ngal = dim(cat.flux)[0];

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

    if (n_elements(cat.ra) == 0 || n_elements(cat.dec) == 0) {
        cat.ra = cat.pos.ra;
        cat.dec = cat.pos.dec;
    }

    if (n_elements(cat.ra) == 0 || n_elements(cat.dec) == 0) {
        print("error: no RA/Dec positions in this file (expected RA, DEC or POS.RA, POS.DEC)");
        return 1;
    }

    print("Photometry: [", ngal, " sources]");
    vec1u gidg = where(finite(cat.flux) && finite(cat.flux_err) && cat.flux > 0 && cat.flux_err > 0
        && cat.flux/cat.flux_err > 3);

    vec1u rndid;
    if (gidg.size() > 2000) {
        auto seed = make_seed(42);
        rndid = shuffle(gidg, seed)[uindgen(2000)];
    } else {
        rndid = uindgen(gidg.size());
    }

    std::string hband = "band";
    std::string hnote = "note";
    std::string hlam  = "lambda [um]";
    std::string hdepuJy = "depth [uJy]";
    std::string hdepAB = "depth [AB]";
    std::string hdet = "detections";
    std::string hd3s = "3-sigma";
    std::string hd5s = "5-sigma";
    std::string har = "area [deg^2]";

    uint_t maxb = std::max(hband.size(), max(length(cat.bands)));
    uint_t maxn = std::max(hnote.size(), max(length(cat.notes)));
    uint_t maxl = std::max(hlam.size(), max(length(strna(cat.lambda))));
    uint_t maxdJ = std::max(hdepuJy.size(), max(length(strna_sci(cat.flux_err[rndid]))));
    uint_t maxdAB = std::max(hdepAB.size(), max(length(strna(cat.flux_err[rndid]))));
    uint_t maxdet = std::max(hdet.size(), length(strn(ngal)));
    uint_t maxd3 = std::max(hd3s.size(), length(strn(ngal)));
    uint_t maxd5 = std::max(hd5s.size(), length(strn(ngal)));
    uint_t maxar = std::max(har.size(), max(length(strna(cat.lambda))));

    std::string header = " "+
        align_center(hband, maxb)+" | "+
        align_center(hnote, maxn)+" | "+
        align_center(hlam, maxl)+" | "+
        align_center(hdepuJy, maxdJ)+" | "+
        align_center(hdepAB, maxdAB)+" | "+
        align_center(hdet, maxdet)+" | "+
        align_center(hd3s, maxd3)+" | "+
        align_center(hd5s, maxd5)+" | "+
        align_center(har, maxar)+" ";

    print(header);
    print(std::string(header.size(), '='));

    for (uint_t i = 0; i < n_elements(cat.bands); ++i) {
        auto f = cat.flux(_,i);
        auto e = cat.flux_err(_,i);

        vec1u idg = where(finite(f) && finite(e) && f > 0 && e > 0);
        vec1u idg3s = where(finite(f) && finite(e) && f > 0 && e > 0 && f/e > 3);

        if (idg.empty() || idg3s.empty()) continue;

        print(" ",
            align_center(cat.bands[i], maxb), " | ",
            align_center(cat.notes[i], maxn), " | ",
            align_right(strn(cat.lambda[i]), maxl), " | ",
            align_right(strn_sci(3*median(e[idg3s])), maxdJ), " | ",
            align_right(strn(uJy2mag(3*median(e[idg3s]))), maxdAB), " | ",
            align_right(strn(idg.size()), maxdet), " | ",
            align_right(strn(total(f[idg]/e[idg] > 3)), maxd3), " | ",
            align_right(strn(total(f[idg]/e[idg] > 5)), maxd5), " | ",
            align_right(strn(field_area(cat.ra[idg3s], cat.dec[idg3s])), maxar)
        );
    }

    return 0;
}
