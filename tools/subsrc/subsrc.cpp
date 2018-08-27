#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string bands_notes, bands, notes;
    std::string psf_model;
    vec1u include, exclude;
    std::string out, reg;
    double fconv = 1.0;
    read_args(argc-2, argv+2, arg_list(
        bands_notes, bands, notes, name(psf_model, "psf"), include, exclude, out, fconv, reg
    ));

    if (psf_model.empty()) {
        error("missing PSF model: psf=...");
        return 1;
    }

    fits::header hdr;
    vec2d img;
    fits::read(argv[1], img, hdr);

    struct {
        vec1d ra, dec;
        vec1f flux;
    } cat;

    std::string cat_file = argv[2];
    if (ends_with(cat_file, ".fits")) {
        struct {
            vec1u id;
            vec1d ra, dec;
            vec1s bands, notes;
        } fcat;

        fits::read_table_loose(cat_file, ftable(fcat.id, fcat.ra, fcat.dec, fcat.bands, fcat.notes));
        if (fcat.ra.empty()) {
            fits::read_table_loose(cat_file, "pos.ra", fcat.ra, "pos.dec", fcat.dec);
            if (fcat.ra.empty()) {
                error("could not find position in this catalog (either RA/DEC or POS.RA/POS.DEC)");
                return 1;
            }
        }

        if (fcat.id.empty()) {
            fcat.id = indgen(fcat.ra.size());
        }

        vec1f flux;
        bool found = false;
        auto info = fits::read_table_columns(cat_file);
        for (auto& i : info) {
            i.name = to_lower(i.name);
            if (i.name == "flux") {
                if (i.dims.size() == 2) {
                    vec2f tflux;
                    fits::read_table(cat_file, "flux", tflux);
                    vec1b bm, nm;

                    if (bands_notes.empty()) {
                        if (bands.empty() || fcat.bands.empty()) {
                            bm = replicate(true, tflux.dims[1]);
                        } else {
                            bm = regex_match(fcat.bands, bands);
                        }
                        if (notes.empty() || fcat.notes.empty()) {
                            nm = replicate(true, tflux.dims[1]);
                        } else {
                            nm = regex_match(fcat.notes, notes);
                        }
                    } else {
                        if (fcat.bands.empty()) {
                            bm = replicate(true, tflux.dims[1]);
                        } else if (fcat.notes.empty()) {
                            bm = regex_match(fcat.bands, bands_notes);
                        } else {
                            bm = regex_match(fcat.bands+" "+fcat.notes, bands_notes);
                        }

                        nm = replicate(true, tflux.dims[1]);
                    }

                    vec1u ib = where(nm && bm);
                    if (ib.empty()) {
                        error("could not find band");
                        return 1;
                    } else if (ib.size() != 1) {
                        error("multiple matches for band");
                        note(fcat.bands[ib]+" "+fcat.notes[ib]);
                        return 1;
                    }

                    flux = tflux(_,ib[0]);
                } else if (i.dims.size() == 1) {
                    fits::read_table(cat_file, "flux", flux);
                } else {
                    error("FLUX column must be 1 or 2 dimensional (got: ", i.dims, ")");
                    return 1;
                }
                found = true;
                break;
            }
        }

        if (!found) {
            error("could not find flux in this catalog (FLUX column)");
            return 1;
        }

        vec1u id = where(
            (include.empty() ? !is_any_of(fcat.id, exclude) : is_any_of(fcat.id, include)) &&
            is_finite(flux) && flux > 0.0
        );

        cat.ra = fcat.ra[id];
        cat.dec = fcat.dec[id];
        cat.flux = flux[id];
    } else {
        vec1f flx;
        ascii::read_table(cat_file, cat.ra, cat.dec, cat.flux);
    }

    print(cat.ra.size(), " sources to subtract from the image");

    astro::wcs w(hdr);
    vec1d x, y;
    astro::ad2xy(w, cat.ra, cat.dec, x, y);

    vec2d psf(img.dims);
    auto pg = progress_start(cat.ra.size());
    for (uint_t i : range(cat.ra)) {
        if (!make_psf({{img.dims[0], img.dims[1]}}, y[i]-1, x[i]-1, psf_model, psf)) {
            return 1;
        }

        img -= psf*cat.flux[i]/fconv;
        progress(pg);
    }

    if (out.empty()) {
        out = file::remove_extension(argv[1])+"_res.fits";
    }

    fits::write(out, img, hdr);

    if (!reg.empty()) {
        print_radec(reg, cat.ra, cat.dec);
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("subsrc v1.0");
    paragraph("usage: subsrc img.fits srcs.[fits/cat] [psf=...,out=...]");
}
