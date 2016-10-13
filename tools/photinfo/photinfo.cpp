#include <phypp.hpp>
#include <phypp/astro/qxmatch.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc == 1) {
        return 0;
    }

    read_args(argc-1, argv+1, arg_list(data_dir));

    auto fdb = read_filter_db(data_dir+"/fits/filters/db.dat");

    if (std::string(argv[1]) == "list") {
        print_filters(fdb);
        return 0;
    }

    struct {
        vec1d ra, dec;
        vec2f flux, flux_err;
        vec1s bands, notes;
        vec1f lambda;
    } cat;

    fits::read_table_loose(argv[1], ftable(
        cat.ra, cat.dec, cat.flux, cat.flux_err, cat.bands, cat.notes, cat.lambda
    ));

    uint_t ngal = cat.flux.dims[0];

    if (cat.lambda.empty()) {
        cat.lambda = fltarr(cat.bands.size());
        for (uint_t i = 0; i < cat.bands.size(); ++i) {
            filter_t filter;
            if (get_filter(fdb, cat.bands[i], filter)) {
                cat.lambda[i] = filter.rlam;
            } else {
                cat.lambda[i] = fnan;
            }
        }
    }

    if (cat.notes.empty()) {
        cat.notes = replicate("", cat.bands.size());
    }

    if (cat.ra.empty() || cat.dec.empty()) {
        fits::read_table_loose(argv[1], "pos.ra", cat.ra, "pos.dec", cat.dec);
    }

    if (cat.ra.empty() || cat.dec.empty()) {
        print("error: no RA/Dec positions in this file (expected RA, DEC or POS.RA, POS.DEC)");
        return 1;
    }

    print("Photometry: [", ngal, " sources]");
    vec1u gidg = where(is_finite(cat.flux) && is_finite(cat.flux_err) && cat.flux > 0 && cat.flux_err > 0
        && cat.flux/cat.flux_err > 3);

    auto seed = make_seed(42);

    vec1u rndid;
    if (gidg.size() >= 2000) {
        rndid = shuffle(seed, gidg)[_-2000];
    } else {
        rndid = uindgen(gidg.size());
    }

    vec1f areas(cat.bands.size());
    vec1f mindist(cat.bands.size());

    std::string cache = argv[1];
    if (end_with(cache, ".fits")) cache = erase_end(cache, 5);
    cache = "."+cache+"_photinfo_cache.fits";

    if (!file::exists(cache) || file::is_older(cache, argv[1])) {
        auto pg = progress_start(cat.bands.size());
        for (uint_t i = 0; i < cat.bands.size(); ++i) {
            auto f = cat.flux(_,i);
            auto e = cat.flux_err(_,i);

            vec1u idg3s = where(is_finite(f) && is_finite(e) && f > 0 && e > 0 && f/e > 3);
            if (!idg3s.empty()) {
                auto hull = build_convex_hull(cat.ra[idg3s], cat.dec[idg3s]);
                areas[i] = field_area_hull(hull);

                if (idg3s.size() > 1000) {
                    double mx = 0.5*(max(hull.x) + min(hull.x));
                    double my = 0.5*(max(hull.y) + min(hull.y));
                    hull.x = (hull.x - mx)*sqrt(1000.0/idg3s.size()) + mx;
                    hull.y = (hull.y - my)*sqrt(1000.0/idg3s.size()) + my;
                }

                vec1u central = idg3s[where(in_convex_hull(cat.ra[idg3s], cat.dec[idg3s], hull))];
                if (!central.empty()) {
                    qxmatch_params p; p.thread = central.size() > 800 ? 2 : 1;
                    auto res = qxmatch(cat.ra[central], cat.dec[central], p);
                    mindist[i] = min(res.d);
                }
            }

            progress(pg);
        }

        fits::write_table(cache, ftable(areas, mindist));
    } else {
        fits::read_table(cache, ftable(areas, mindist));
    }

    std::string hid = "id";
    std::string hband = "band";
    std::string hnote = "note";
    std::string hlam  = "lambda [um]";
    std::string hdepuJy = "depth [uJy]";
    std::string hdepAB = "depth [AB]";
    std::string hdet = "detections";
    std::string hd3s = "3-sigma";
    std::string hd5s = "5-sigma";
    std::string har = "area [deg^2]";
    std::string hmd = "min dist [\"]";
    std::string had = "avg dist [\"]";

    uint_t maxi = std::max(hid.size(), uint_t(1+log10(cat.bands.size())));
    uint_t maxb = std::max(hband.size(), max(length(cat.bands)));
    uint_t maxn = std::max(hnote.size(), max(length(cat.notes)));
    uint_t maxl = std::max(hlam.size(), max(length(strna(cat.lambda))));
    uint_t maxdJ = std::max(hdepuJy.size(), max(length(strna_sci(cat.flux_err[rndid]))));
    uint_t maxdAB = std::max(hdepAB.size(), max(length(strna(cat.flux_err[rndid]))));
    uint_t maxdet = std::max(hdet.size(), length(strn(ngal)));
    uint_t maxd3 = std::max(hd3s.size(), length(strn(ngal)));
    uint_t maxd5 = std::max(hd5s.size(), length(strn(ngal)));
    uint_t maxar = std::max(har.size(), max(length(strna(cat.lambda))));
    uint_t maxmd = std::max(hmd.size(), max(length(strna(mindist))));
    uint_t maxad = std::max(had.size(), max(length(strna(mindist))));

    std::string header = " "+
        align_center(hid, maxi)+" | "+
        align_center(hband, maxb)+" | "+
        align_center(hnote, maxn)+" | "+
        align_center(hlam, maxl)+" | "+
        align_center(hdepuJy, maxdJ)+" | "+
        align_center(hdepAB, maxdAB)+" | "+
        align_center(hdet, maxdet)+" | "+
        align_center(hd3s, maxd3)+" | "+
        align_center(hd5s, maxd5)+" | "+
        align_center(har, maxar)+" | "+
        align_center(hmd, maxmd)+" | "+
        align_center(had, maxad)+" ";

    print(header);
    print(std::string(header.size(), '='));

    for (uint_t i = 0; i < cat.bands.size(); ++i) {
        auto f = cat.flux(_,i);
        auto e = cat.flux_err(_,i);

        vec1u idg = where(is_finite(f) && is_finite(e) && f > 0 && e > 0);
        vec1u idg3s = where(is_finite(f) && is_finite(e) && f > 0 && e > 0 && f/e > 3);

        if (idg.empty() || idg3s.empty()) continue;

        print(" ",
            align_center(strn(i), maxi), " | ",
            align_center(cat.bands[i], maxb), " | ",
            align_center(cat.notes[i], maxn), " | ",
            align_right(strn(cat.lambda[i]), maxl), " | ",
            align_right(strn_sci(3*median(e[idg3s])), maxdJ), " | ",
            align_right(keep_start(strn(float(uJy2mag(3*median(e[idg3s])))), 4), maxdAB), " | ",
            align_right(strn(idg.size()), maxdet), " | ",
            align_right(strn(total(f[idg]/e[idg] > 3)), maxd3), " | ",
            align_right(strn(total(f[idg]/e[idg] > 5)), maxd5), " | ",
            align_right(align_left(strn(areas[i]), 9, '0'), maxar), " | ",
            align_right(strn(mindist[i]), maxmd), " | ",
            align_right(strn(3600.0*sqrt(areas[i]/idg3s.size())), maxad)
        );
    }

    return 0;
}
