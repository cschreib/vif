#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string bands_notes, bands, notes;
    std::string psf_model;
    vec1u include, exclude;
    std::string out, reg;
    double map_unit = 1.0;
    read_args(argc-2, argv+2, arg_list(
        bands_notes, bands, notes, name(psf_model, "psf"), include, exclude, out, map_unit, reg
    ));

    if (psf_model.empty()) {
        error("missing PSF model");
        return 1;
    }

    fits::header hdr;
    vec2d img;
    fits::read(argv[1], img, hdr);

    struct {
        vec1u id;
        vec1d ra, dec;
        vec2f flux;
        vec1s bands, notes;
    } cat;

    fits::read_table_loose(argv[2], cat);
    if (cat.ra.empty()) {
        fits::read_table_loose(argv[2], "pos.ra", cat.ra, "pos.dec", cat.dec);
        if (cat.ra.empty()) {
            error("could not find position in this catalog (either RA/DEC or POS.RA/POS.DEC)");
            return 1;
        }
    }

    if (cat.id.empty()) {
        cat.id = uindgen(cat.ra.size());
    }

    vec1b bm, nm;
    if (bands_notes.empty()) {
        if (bands.empty() || cat.bands.empty()) {
            bm = replicate(true, cat.flux.dims[1]);
        } else {
            bm = match(cat.bands, bands);
        }
        if (notes.empty() || cat.notes.empty()) {
            nm = replicate(true, cat.flux.dims[1]);
        } else {
            nm = match(cat.notes, notes);
        }
    } else {
        if (cat.bands.empty()) {
            bm = replicate(true, cat.flux.dims[1]);
        } else if (cat.notes.empty()) {
            bm = match(cat.bands, bands_notes);
        } else {
            bm = match(cat.bands+" "+cat.notes, bands_notes);
        }

        nm = replicate(true, cat.flux.dims[1]);
    }

    vec1u ib = where(nm && bm);
    if (ib.empty()) {
        error("could not find band");
        return 1;
    } else if (ib.size() != 1) {
        error("multiple matches for band");
        note(strn(cat.bands[ib]+" "+cat.notes[ib]));
        return 1;
    }

    vec1u id = where(
        (include.empty() ? !is_any_of(cat.id, exclude) : is_any_of(cat.id, include)) &&
        finite(cat.flux(_,ib[0])) && cat.flux(_,ib[0]) > 0.0
    );

    print(id.size(), " sources to subtract from the image");

    fits::wcs astro(hdr);
    vec1d x, y;
    fits::ad2xy(astro, cat.ra[id], cat.dec[id], x, y);

    vec2d psf(img.dims);
    auto pg = progress_start(id.size());
    for (uint_t i : range(id)) {
        if (!make_psf({{img.dims[0], img.dims[1]}}, y[i]-1, x[i]-1, psf_model, psf)) {
            return 1;
        }

        img -= map_unit*psf*cat.flux(id[i],ib[0]);
        progress(pg);
    }

    if (out.empty()) {
        out = erase_end(argv[1], 5)+"_res.fits";
    }

    fits::write(out, img, hdr);

    if (!reg.empty()) {
        print_radec(reg, cat.ra[id], cat.dec[id]);
    }

    return 0;
}

void print_help() {
    print("[WIP]");
}
